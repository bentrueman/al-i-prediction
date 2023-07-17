
source("R/ridge-t/model-combined.R")
library("ggdist")

# ------------------------ impute  ------------------------

these_cens <- get_cens(sdata, these_vars_noy)

data_test_imputed_cens <- impute(
  data_combined, model, 
  var = these_vars_noy,
  mi = NULL, 
  cens = these_cens
) %>% 
  # remove training data:
  map(~ filter(.x, type == "test_cens")) %>% 
  # add uncensored rows of test data to each list element:
  map(~ bind_rows(.x, filter(data_test, !rowid %in% data_test_subset$rowid))) %>% 
  # arrange by rowid (to match data_test):
  map(~ arrange(.x, rowid))

# ------------------------ predict  ------------------------

# go back to trained model:

source("R/ridge-t/model-train.R")

# inputs:

post_draws <- as_draws_df(model)

r_site <- post_draws %>% 
  select(starts_with("r_1_aliugl_1")) %>% 
  t()

gam_preds <- posterior_smooths(
  model_brms, 
  newdata = data_test_imputed_cens[[1]], 
  smooth = "s(yday, bs = \"cc\")",
  resp = "tempdegcsonde"
)

# make site index for test data:

site_id <- tibble(
  site = data_train$site, 
  id = as.numeric(sdata$J_1_aliugl)
) %>% 
  distinct(site, id) %>% 
  right_join(data_test_imputed_cens[[1]], by = "site") %>% 
  arrange(rowid) %>% 
  select(rowid, id) %>% 
  verify(rowid == seq_len(nrow(data_test)))

# impute missings via posterior prediction:

data_test_imputed <- map2(
  data_test_imputed_cens,
  seq_len(nrow(post_draws)),
  ~ .x %>% 
    mutate(
      # first, impute variables that don't depend on any other variables:
      # DOC:
      avg_doc_mg_l = if_else(
        is.na(avg_doc_mg_l), 
        post_draws$Intercept_avgdocmgl[.y],
        avg_doc_mg_l
      ),
      # Ca:
      calcium_dis_mg_l = if_else(
        is.na(calcium_dis_mg_l), 
        post_draws$Intercept_calciumdismgl[.y],
        calcium_dis_mg_l
      ),
      # fluoride:
      fluoride_ug_l = if_else(
        is.na(fluoride_ug_l), 
        post_draws$Intercept_fluorideugl[.y],
        fluoride_ug_l
      ),
      # dissolved aluminum:
      aluminum_dis_ug_l = if_else(
        is.na(aluminum_dis_ug_l), 
        post_draws$`bsp_aluminumdisugl[1]`[.y] * avg_doc_mg_l,
        aluminum_dis_ug_l
      ),
      # dissolved titanium:
      titanium_dis_ug_l = if_else(
        is.na(titanium_dis_ug_l), 
        post_draws$`bsp_titaniumdisugl[1]`[.y] * aluminum_dis_ug_l,
        titanium_dis_ug_l
      ),
      # colour:
      avg_colour_colour_units = if_else(
        is.na(avg_colour_colour_units), 
        post_draws$`bsp_avgcolourcolourunits[1]`[.y] * avg_doc_mg_l,
        avg_colour_colour_units
      ),
      # dissolved iron:
      iron_dis_ug_l = if_else(
        is.na(iron_dis_ug_l), 
        post_draws$`bsp_irondisugl[1]`[.y] * avg_doc_mg_l,
        iron_dis_ug_l
      ),
      # dissolved cerium:
      cerium_dis_ug_l = if_else(
        is.na(cerium_dis_ug_l), 
        post_draws$`bsp_ceriumdisugl[1]`[.y] * iron_dis_ug_l,
        cerium_dis_ug_l
      ),
      # pH (sonde) :
      ph_sonde = if_else(
        is.na(ph_sonde), 
        post_draws$`bsp_phsonde[1]`[.y] * aluminum_dis_ug_l,
        ph_sonde
      ),
      # water temperature (sonde):
      temp_deg_c_sonde  = if_else(
        is.na(temp_deg_c_sonde), 
        gam_preds[.y, ], 
        temp_deg_c_sonde
      ),
      # sulfate
      sulfate_mg_l = if_else(
        is.na(sulfate_mg_l), 
        post_draws$`bsp_sulfatemgl[1]`[.y] * calcium_dis_mg_l,
        sulfate_mg_l
      ),
      # alkalinity
      alkalinity_mg_ca_co3l = if_else(
        is.na(alkalinity_mg_ca_co3l), 
        post_draws$`bsp_alkalinitymgcaco3l[1]`[.y] * calcium_dis_mg_l,
        alkalinity_mg_ca_co3l
      ),
      # lithium
      lithium_dis_ug_l = if_else(
        is.na(lithium_dis_ug_l), 
        post_draws$`bsp_lithiumdisugl[1]`[.y] * fluoride_ug_l,
        lithium_dis_ug_l
      )
    ), .progress = TRUE 
  )

# ------------------------ predict  ------------------------

yhat <- predict_stan(
  post_draws, r_site, data_test_imputed, data_scaled, data_test_raw,
  site_id$id
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

write_csv(yhat$yhat_summ, "data-clean/yhat-test.csv")

bind_rows(
  rmse = errors$rmse,
  mae = errors$mae,
  .id = "metric"
) %>% 
  write_csv("data-clean/error-test.csv")
