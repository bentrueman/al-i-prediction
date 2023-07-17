
source("R/ridge-t/model-setup.R")

# ------------------------ read  ------------------------

# exclude response:
these_vars_noy <- these_vars$var[these_vars$var != "ali_ug_l"] 
cens_ind_noy <- cens_ind[cens_ind != "cens_ali_ug_l"]

# combine the training set with those rows of the test set that have censored values:

data_test_raw <- read_csv("data-clean/data-al-i-test.csv")

# scale test data using col means and sds from training data:

trainscale <- list(
  mean = attr(data_scaled, "scaled:center")[-1], 
  sd = attr(data_scaled, "scaled:scale")[-1]
) %>% 
  map(
    ~ matrix(
      .x, 
      nrow = nrow(data_test_raw), 
      ncol = length(these_vars_noy),
      byrow = TRUE
    )
  )

data_test_scaled <- (select(data_test_raw, all_of(these_vars_noy)) - trainscale$mean) / 
  trainscale$sd

# censoring indicators:

data_test_cens <- data_test_raw %>% 
  select(all_of(c(cens_ind_noy, "site", "yday")))

# combine:

data_test <- cbind(data_test_scaled, data_test_cens) %>%
  rowid_to_column()

data_test_subset <- data_test %>% 
  # retain only rows of test data with censored values:
  filter(if_any(starts_with("cens_"), ~ .x == "left"))

data_combined <- bind_rows(
  "train" = data_train, 
  "test_cens" = data_test_subset, 
  .id = "type"
)
  
# ------------------------ stan inputs ------------------------

# id variables with censoring:

pcens <- id_cens(data_combined)

sdata <- make_standata(bform, data = data_combined, knots = list(yday = 0:1))

for(i in 1:length(pcens$name)) {
  var <- pcens$name[i]
  sdata <- modify_standata(sdata, data_combined, dls[[var]], var)
}

scode <- make_stancode(
  bform, 
  data = data_combined, 
  knots = list(yday = 0:1), 
  stanvars = stanvars,
  family = "student",
  prior = bprior
) %>% 
  modify_stancode(pcens$name)

# ------------------------ fit stan model ------------------------

bname <- "model-t-combined"

model <- fit_model("models", bname, scode, sdata, stanseed, max_treedepth = 14)
