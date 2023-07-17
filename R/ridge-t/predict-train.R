
source("R/ridge-t/model-train.R")
library("ggdist")
library("furrr")
plan("multisession")

# ------------------------ impute  ------------------------

these_vars_noy <- these_vars$var[these_vars$var != "ali_ug_l"] # exclude response

these_mi <- get_mi(sdata, these_vars_noy)

these_cens <- get_cens(sdata, these_vars_noy)

data_train_imputed <- impute(
  data_train, 
  model = model, 
  var = these_vars_noy,
  mi = these_mi, 
  cens = these_cens
)

# ------------------------ predict  ------------------------

# group-level intercepts:

post_draws <- as_draws_df(model)

r_site <- post_draws %>% 
  select(starts_with("r_1_aliugl_1")) %>% 
  t()

# predict Ali:

yhat <- predict_stan(
  post_draws, r_site, data_train_imputed, data_scaled, data_train_raw,
  sdata$J_1_aliugl
)

# ------------------------ visualize  ------------------------

yhat$yhat_summ %>% 
  ggplot(aes(ali_ug_l, y)) + 
  geom_abline(col = "red") +
  geom_point(
    data = . %>% 
      filter(cens_ali_ug_l != "left")
  ) + 
  geom_segment(
    data = . %>% 
      filter(cens_ali_ug_l == "left"),
    aes(xend = -Inf, yend = y)
  )

# ------------------------ error  ------------------------

errors <- estimate_error(yhat)

# ------------------------ write  ------------------------

write_csv(yhat$yhat_summ, "data-clean/yhat-train.csv")

bind_rows(
  rmse = errors$rmse,
  mae = errors$mae,
  .id = "metric"
) %>% 
  write_csv("data-clean/error-train.csv")
