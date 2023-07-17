
library("tidyverse")
library("janitor", include.only = "clean_names")
library("brms")
library("cmdstanr")
library("assertr")
source("R/ridge-t/functions.R")

# ------------------------ select vars ------------------------

these_vars <- tibble::tribble(
                       ~var,               ~varstan, ~predict,
                 "ali_ug_l",               "aliugl",    "ali",
        "aluminum_dis_ug_l",       "aluminumdisugl",    "ali",
        "titanium_dis_ug_l",       "titaniumdisugl",    "ali",
             "avg_doc_mg_l",            "avgdocmgl",    "ali",
  "avg_colour_colour_units", "avgcolourcolourunits",    "ali",
            "iron_dis_ug_l",           "irondisugl",    "ali",
          "cerium_dis_ug_l",         "ceriumdisugl",    "ali",
                 "ph_sonde",              "phsonde",    "ali",
         "temp_deg_c_sonde",        "tempdegcsonde",    "ali",
             "sulfate_mg_l",           "sulfatemgl",    "ali",
    "alkalinity_mg_ca_co3l",   "alkalinitymgcaco3l",    "ali",
         "lithium_dis_ug_l",        "lithiumdisugl",    "ali",
         "calcium_dis_mg_l",        "calciumdismgl",     "mi",
            "fluoride_ug_l",          "fluorideugl",     "mi"
  )

cens_ind <- paste0("cens_", these_vars$var)

# ------------------------ prep data ------------------------

data_train_raw <- read_csv("data-clean/data-al-i-train.csv")

# censoring indicators:

data_cens <- data_train_raw %>% 
  select(all_of(c(cens_ind, "site", "yday")))

# scale predictors:

data_scaled <- data_train_raw %>% 
  select(all_of(these_vars$var)) %>% 
  scale()

# combine:

data_train <- cbind(data_scaled, data_cens)

# detection limits:

dls <- read_csv("data-clean/detection-limits-modified.csv") %>% 
  pivot_wider(names_from = parameter, values_from = detection_limit) %>% 
  clean_names()

# id variables with censoring:

pcens <- id_cens(data_train)

# ------------------------ options  ------------------------

stanseed <- 1257
options(mc.cores = parallel::detectCores())
  
# ------------------------ formula  ------------------------

bform <- bf(ali_ug_l | mi() ~ 0 +
              mi(aluminum_dis_ug_l) + 
              mi(titanium_dis_ug_l) +
              mi(avg_doc_mg_l) + 
              mi(avg_colour_colour_units) +
              mi(iron_dis_ug_l) +
              mi(cerium_dis_ug_l) +
              mi(ph_sonde) + 
              mi(temp_deg_c_sonde) +
              mi(sulfate_mg_l) +
              mi(alkalinity_mg_ca_co3l) +
              mi(lithium_dis_ug_l) +
              mi(calcium_dis_mg_l) +
              mi(fluoride_ug_l) +
              (1 | site)) + 
  # for each missing value model, pick the best predictor based on 
  # Pearson correlation in the training set.
  # the exceptions are Ca, DOC, and F which are predicted 
  # using intercept-only models.
  bf(aluminum_dis_ug_l | mi() ~ 0 + mi(avg_doc_mg_l)) +
  bf(titanium_dis_ug_l | mi() ~ 0 + mi(aluminum_dis_ug_l)) +
  bf(avg_colour_colour_units | mi() ~ 0 + mi(avg_doc_mg_l)) +
  bf(avg_doc_mg_l | mi() ~ 1) +
  bf(iron_dis_ug_l | mi() ~ 0 + mi(avg_doc_mg_l)) +
  bf(cerium_dis_ug_l | mi() ~ 0 + mi(iron_dis_ug_l)) +
  bf(ph_sonde | mi() ~ 0 + mi(aluminum_dis_ug_l)) + 
  bf(calcium_dis_mg_l | mi() ~ 1) +
  bf(temp_deg_c_sonde | mi() ~ 0 + s(yday, bs = "cc")) +
  bf(sulfate_mg_l | mi() ~ 0 + mi(calcium_dis_mg_l)) +
  bf(alkalinity_mg_ca_co3l | mi() ~ 0 + mi(calcium_dis_mg_l)) +
  bf(lithium_dis_ug_l | mi() ~ 0 + mi(fluoride_ug_l)) +
  bf(fluoride_ug_l | mi() ~ 1) +
  set_rescor(FALSE)

# ------------------------ stan inputs ------------------------

sdata <- make_standata(bform, data = data_train, knots = list(yday = 0:1))

# add censored observations as parameters:

for(i in 1:length(pcens$name)) {
  var <- pcens$name[i]
  sdata <- modify_standata(sdata, data_train, dls[[var]], var)
}

# set priors:

bprior <- set_prior("normal(0, tau)", class = "b", resp = these_vars$varstan[1]) +
  # 4, 9, 13, and 14 don't have class "b" coefficients
  # ("avgdocmgl", "tempdegcsonde", "calciumdismgl", "fluorideugl", 
  # i.e., these_vars$varstan[c(4,9, 13, 14)])
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[2]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[3]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[5]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[6]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[7]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[8]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[10]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[11]) +
  set_prior("normal(0, 1)", class = "b", resp = these_vars$varstan[12]) +
  # all variables are centered, so intercepts get a tight prior around 0:
  set_prior("normal(0, 0.1)", class = "Intercept", resp = these_vars$varstan[4]) +
  set_prior("normal(0, 0.1)", class = "Intercept", resp = these_vars$varstan[13]) +
  set_prior("normal(0, 0.1)", class = "Intercept", resp = these_vars$varstan[14]) +
  set_prior("target += cauchy_lpdf(tau | 0, 1) - cauchy_lccdf(0 | 0, 1)", check = FALSE)

# add ridge scale parameter:

stanvars <- stanvar(scode = "real<lower=0> tau;", block = "parameters")

# generate stan code:

scode <- make_stancode(
  bform, 
  data = data_train, 
  knots = list(yday = 0:1), 
  stanvars = stanvars,
  family = "student",
  prior = bprior
) %>% 
  modify_stancode(pcens$name)
