
source("R/ridge-t/model-setup.R")

# ------------------------ fit stan model ------------------------

path <- "models"
bname <- "model-t-train"
file <- paste(path, bname, sep = "/")

model <- fit_model(
  path, bname, scode, sdata, stanseed, max_treedepth = 14
)

# store model as a brmsfit object for later prediction of spline terms:

model_rds <- list.files(path, pattern = paste0(bname, "\\.rds"))

if (length(model_rds) > 0) {
  model_brms <- brm(
    formula = bform, 
    data = data_train, 
    knots = list(yday = 0:1), 
    stanvars = stanvars,
    family = "student",
    prior = bprior,
    file = file
  )
} else {
  model_rstan <- fit_model(
    path, bname, scode, sdata, stanseed, rstan = TRUE
  )
  model_brms <- brm(
    formula = bform, 
    data = data_train, 
    knots = list(yday = 0:1), 
    stanvars = stanvars,
    family = "student",
    prior = bprior,
    empty = TRUE
  )
  model_brms$fit <- model_rstan
  model_brms <- rename_pars(model_brms)
  brms:::write_brmsfit(model_brms, file = file)
}

