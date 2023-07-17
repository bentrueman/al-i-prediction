
library("glue")
library("testthat")
library("rstan", include.only = "read_stan_csv")

# ------------------------ modify_standata() ------------------------

modify_standata <- function(sdata, data, lcl, var) {
  
  if (length(lcl) != length(var)) stop("lengths of 'var' and 'lcl' must be equal")
  
  varstan <- str_remove_all(var, "_")
  
  for(i in seq_len(length(var))) {
    # make logical censoring indicator:
    cens_logical <- data[[paste0("cens_", var[i])]] == "left"
    
    if (sum(cens_logical) > 0) { # append to standata list only if there are left-censored values
      sdata[[paste0("Ncens_", varstan[i])]] <- sum(cens_logical) # number of left-censored
      # positions of left-censored:
      sdata[[paste0("Jcens_", varstan[i])]] <- as.array(seq_len(nrow(data))[cens_logical]) 
      sdata[[paste0("U_", varstan[i])]] <- lcl[i] # left-censoring limit
    } else message(glue("No left-censored {var} values."))
  }
  
  sdata
}

test_that("modify_standata() creates the expected list", {
  data <- tibble(x = rep(NA, 11), cens_x = c(rep("none", 10), "left"))
  out <- modify_standata(sdata = list(), data, 999, "x")
  expect_equal(list(Ncens_x = 1, Jcens_x = as.array(11), U_x = 999), out)
})

# ------------------------ modify_stancode() ------------------------

modify_stancode <- function(scode, var) {
  
  var <- str_remove_all(var, "_")
  
  for(i in seq_len(length(var))) {
    
    # modifications to data block:
    n_cens <- glue("int<lower=0> Ncens_{var[i]};  // number of left-censored")
    j_cens <- glue("int<lower=1> Jcens_{var[i]}[Ncens_{var[i]}];  // positions of left-censored")
    u <- glue("real U_{var[i]};  // left-censoring limit")
    # modifications to parameters block:
    y_cens <- glue("vector<upper=U_{var[i]}>[Ncens_{var[i]}] Ycens_{var[i]};  // estimated left-censored")
    # modifications to model block:
    yl <- glue("Yl_{var[i]}[Jcens_{var[i]}] = Ycens_{var[i]}; // add imputed left-censored values")
    
    scode <- scode %>%
      # modifications to data block:
      str_replace(
        "(data \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1", n_cens, j_cens, u), collapse = "\n  ")
      ) %>% 
      # modifications to parameters block:
      str_replace(
        "(parameters \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1\n  ", y_cens), collapse = "")
      ) %>% 
      # modifications to model block:
      str_replace(
        "(model \\{\n(.|\n)*?)(?=\n    mu_)",
        paste(c("\\1\n    ", yl), collapse = "")
      )
    
  }
  
  class(scode) <- "brmsmodel"
  
  scode
  
}

# ------------------------ impute() ------------------------

impute_censored <- function(cens, post_draws, varstan, id, data_imputed, drawids, var) {
  
  # impute censored:
  if (!is.null(cens)) {
    # select censored values:
    censored <- post_draws %>% 
      select(starts_with(paste0("Ycens_", varstan))) %>% 
      t()
    # optionally, subset censored:
    if (!is.null(id)) {
      censored <- censored[id, ]
    }
    # impute:
    data_imputed <- map2(
      data_imputed,
      drawids,
      \(x, y) {
        x[cens, var] <- censored[, y]
        return(x)
      }
    )
  }
  
  return(data_imputed)
  
}

test_that("impute_censored() returns expected list", {
  cens <- c(2, 4, 6, 8, 10)
  post_draws <- tibble(Ycens_x = 1:10)
  data_imputed <- map(1:10, ~ tibble(x = rep(NA, 10)))
  data_imputed <- impute_censored(
    cens, post_draws, varstan = "x", id = NULL, 
    data_imputed, drawids = 1:10, var = "x"
  )
  expect_equal(
    unlist(map(data_imputed, ~ unique(na.omit(.x$x)))), 
    post_draws$Ycens_x
  )
})

impute_missing <- function(mi, post_draws, varstan, data_imputed, drawids, var) {
  # impute missing:
  if (!is.null(mi)) {
    # select missing values:
    missing <- post_draws %>% 
      select(starts_with(paste0("Ymi_", varstan))) %>% 
      t()
    data_imputed <- map2(
      data_imputed,
      drawids,
      \(x, y) {
        x[mi, var] <- missing[, y]
        return(x)
      }
    )
  }
  
  return(data_imputed)
  
}

test_that("impute_missing() returns expected list", {
  mi <- c(2, 4, 6, 8, 10)
  post_draws <- tibble(Ymi_x = 1:10)
  data_imputed <- map(1:10, ~ tibble(x = rep(NA, 10)))
  data_imputed <- impute_missing(mi, post_draws, varstan = "x", data_imputed, drawids = 1:10, var = "x")
  expect_equal(
    unlist(map(data_imputed, ~ unique(na.omit(.x$x)))), 
    post_draws$Ymi_x
  )
})

test_that("impute_missing() and impute_censored() work together", {
  cens <- c(2, 4, 6, 8, 10) - 1
  mi <- c(2, 4, 6, 8, 10)
  post_draws <- tibble(Ymi_x = 1:10, Ycens_x = 11:20, Ymi_y = 11:20, Ycens_y = 1:10)
  data_imputed <- map(1:10, ~ tibble(x = rep(NA, 10), y = rep(NA, 10)))
  # impute x:
  data_imputed <- impute_missing(mi, post_draws, varstan = "x", data_imputed, drawids = 1:10, var = "x")
  data_imputed <- impute_censored(
    cens, post_draws, varstan = "x", id = NULL, 
    data_imputed, drawids = 1:10, var = "x"
  )
  expect_equal(
    unlist(map(data_imputed, ~ unique(.x$x[mi]))), 
    post_draws$Ymi_x
  )
  expect_equal(
    unlist(map(data_imputed, ~ unique(.x$x[cens]))), 
    post_draws$Ycens_x
  )
  expect_equal(unlist(map(data_imputed, ~ unique(.x$y))), rep(NA, 10))
  # impute y:
  data_imputed <- impute_missing(mi, post_draws, varstan = "y", data_imputed, drawids = 1:10, var = "y")
  data_imputed <- impute_censored(
    cens, post_draws, varstan = "y", id = NULL, 
    data_imputed, drawids = 1:10, var = "y"
  )
  expect_equal(
    unlist(map(data_imputed, ~ unique(.x$y[mi]))), 
    post_draws$Ymi_y
  )
  expect_equal(
    unlist(map(data_imputed, ~ unique(.x$y[cens]))), 
    post_draws$Ycens_y
  )
})

impute <- function(data, model, var, mi = NULL, cens = NULL, id = NULL, post_draws = NULL) {
  
  varstan <- str_remove_all(var, "_") # variable names in stan
  
  if (is.null(post_draws)) {
    post_draws <- as_draws_df(model) %>% # posterior draws
      as_tibble()
  }
  
  ndraws <- nrow(post_draws) # number of draws
  
  drawids <- seq_len(ndraws) # sequence of draws
  
  data_imputed <- map(drawids, ~ data) # a list of identical datasets to be imputed
  
  # for each variable:
  
  for(i in seq_len(length(var))) {
    
    data_imputed <- impute_censored(cens[[i]], post_draws, varstan[[i]], id, data_imputed, drawids, var[[i]])
    
    data_imputed <- impute_missing(mi[[i]], post_draws, varstan[[i]], data_imputed, drawids, var[[i]])
    
  }
  
  return(data_imputed)
  
}

test_that("impute() returns expected list", {
  pos <- c(2, 4, 6, 8, 10)
  cens <- list(NULL, pos)
  mi <- list(pos - 1, NULL)
  vars <- list("x", "y")
  post_draws <- tibble(Ymi_x = 10:1, Ycens_y = 1:10)
  data <- tibble(x = rep(NA, 10), y = rep(NA, 10))
  out <- impute(data, model = NULL, vars, mi = mi, cens = cens, post_draws = post_draws)
  out_summ <- map(
    out, 
    ~ .x %>% 
      mutate(z = coalesce(x, y)) %>% 
      pull(z) %>% 
      sum()
  ) %>% 
    unlist() %>% 
    unique()
  expect_equal(11 * 5, out_summ)
})

# ------------------------ rmse ------------------------

rmse <- function(y, yhat) {
  sqrt(mean((y - yhat) ^ 2))
}

# ------------------------ mae ------------------------

mae <- function(y, yhat) {
  median(abs(y - yhat))
}

# ------------------------ get_cens() ------------------------

get_cens <- function(sdata, these_vars) {
  these_vars_stan <- str_remove_all(these_vars, "_")
  these_cens <- paste0("Jcens_", str_remove_all(these_vars, "_"))
  sdata[these_cens]
}

# ------------------------ get_mi() ------------------------

get_mi <- function(sdata, these_vars) {
  these_vars_stan <- str_remove_all(these_vars, "_")
  these_mi <- paste0("Jmi_", str_remove_all(these_vars, "_"))
  sdata[these_mi]
}

# ------------------------ id_cens() ------------------------

id_cens <- function(data) {
  data %>% 
    select(starts_with("cens_")) %>% 
    pivot_longer(everything()) %>% 
    group_by(name = str_extract(name, "(?<=cens_).+")) %>% 
    summarize(pcens = mean(value == "left")) %>% 
    ungroup() %>% 
    filter(pcens > 0)
}

# ------------------------ predict_stan() ------------------------

predict_stan <- function(
    post_draws, 
    r_site, 
    data_train_imputed, 
    data_scaled, 
    data_train_raw,
    index) {
  
  draw_ids <- seq_len(nrow(post_draws))
  
  yhat <- map2(
    draw_ids,
    data_train_imputed,
    \(x, y) {
      yhat <- y$aluminum_dis_ug_l * post_draws$`bsp_aliugl[1]`[x] + 
        y$titanium_dis_ug_l * post_draws$`bsp_aliugl[2]`[x] + 
        y$avg_doc_mg_l * post_draws$`bsp_aliugl[3]`[x] + 
        y$avg_colour_colour_units * post_draws$`bsp_aliugl[4]`[x] + 
        y$iron_dis_ug_l * post_draws$`bsp_aliugl[5]`[x] + 
        y$cerium_dis_ug_l * post_draws$`bsp_aliugl[6]`[x] + 
        y$ph_sonde * post_draws$`bsp_aliugl[7]`[x] + 
        y$temp_deg_c_sonde * post_draws$`bsp_aliugl[8]`[x] + 
        y$sulfate_mg_l * post_draws$`bsp_aliugl[9]`[x] + 
        y$alkalinity_mg_ca_co3l * post_draws$`bsp_aliugl[10]`[x] + 
        y$lithium_dis_ug_l * post_draws$`bsp_aliugl[11]`[x] + 
        y$calcium_dis_mg_l * post_draws$`bsp_aliugl[12]`[x] + 
        y$fluoride_ug_l * post_draws$`bsp_aliugl[13]`[x] + 
        r_site[index, x] # site-level intercepts
    }
  )
  
  pp <- map2(
    draw_ids, 
    yhat, 
    \(x, y) {
      rstudent_t(
        y, 
        df = post_draws$nu_aliugl[x], 
        mu = y, 
        sigma = post_draws$sigma_aliugl[x]
      ) * attr(data_scaled, "scaled:scale")[1] + 
        attr(data_scaled, "scaled:center")[1]
    }
  )
  
  # retransform yhat:
  
  yhat <- map(yhat, ~ .x * attr(data_scaled, "scaled:scale")[1] + 
                attr(data_scaled, "scaled:center")[1]) 
  
  yhat_summ <- yhat %>% 
    do.call(rbind, .) %>%
    apply(2, median_qi) %>% 
    bind_rows() %>% 
    as_tibble()
  
  yhat_summ$ali_ug_l <- data_train_raw$ali_ug_l
  yhat_summ$cens_ali_ug_l <- data_train_raw$cens_ali_ug_l
  
  list(yhat = yhat, yhat_summ = yhat_summ, pp = pp)
}

# ------------------------ estimate_error() ------------------------

estimate_error <- function(yhat) {
  
  cens_or_NA <- yhat$yhat_summ$cens_ali_ug_l == "left" | 
    is.na(yhat$yhat_summ$ali_ug_l)
  
  rmse <- map_dbl(
    yhat$yhat, 
    ~ rmse(yhat$yhat_summ$ali_ug_l[!cens_or_NA], pmax(.x[!cens_or_NA], 0))
  )
  
  mae <- map_dbl(
    yhat$yhat,
    ~ mae(yhat$yhat_summ$ali_ug_l[!cens_or_NA], pmax(.x[!cens_or_NA], 0))
  )
  
  list(rmse = median_qi(rmse), mae = median_qi(mae))
  
}

# ------------------------ fit_model() ------------------------

fit_model <- function(path, bname, scode, sdata, stanseed, rstan = FALSE, ...) {
  
  model_csvs <- list.files(
    path, pattern = paste0(bname, "-\\d+-[1-4]\\.csv"), full.names = TRUE
  )
  
  if (length(model_csvs) > 0) {
    if (rstan) {
      # for feeding back into brms:
      model <- read_stan_csv(model_csvs)
    } else {
      # for regular use:
      model <- as_cmdstan_fit(model_csvs) 
    }
  } else {
    # setup model:
    model_setup <- cmdstan_model(stan_file = write_stan_file(scode), compile = FALSE)
    # format:
    model_setup$format(overwrite_file = TRUE, canonicalize = TRUE, backup = FALSE)
    # write formatted code to file:
    write_stan_file(model_setup$code(), dir = "models", basename = bname)
    # compile:
    model_setup$compile()
    # sample:
    model <- model_setup$sample(data = sdata, seed = stanseed, ...)
    # save output:
    model$save_output_files(
      dir = "models", basename = bname, random = FALSE
    )
  }
  
  model
  
}

# ------------------------ predict_pr() ------------------------

predict_pr <- function(
    yhat,
    post_draws,
    data_train_imputed, 
    data_scaled,
    xi,
    beta_i) {
  
  draw_ids <- seq_len(nrow(post_draws)) 
  
  # full model residuals
  
  r <- map(
    yhat,
    # rescale:
    ~ (.x - attr(data_scaled, "scaled:center")[1]) / 
      attr(data_scaled, "scaled:scale")[1] - 
      data_scaled[, "ali_ug_l"]
  ) 
  
  beta_xi <- map2(
    draw_ids,
    data_train_imputed,
    \(x, y) {
      yhat <- y[[xi]] * post_draws[[beta_i]][x]
    }
  )
  
  pr <- map2(r, beta_xi, ~ .x + .y)
  
  # summarize partial residuals:
  
  pr_summ <- pr %>% 
    do.call(rbind, .) %>%
    apply(2, median_qi) %>% 
    bind_rows() %>% 
    as_tibble()
  
  list(pr = pr, pr_summ = pr_summ)
}
