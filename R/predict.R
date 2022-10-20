
#' Estimates treatment effects by regressing outcome residuals on treatment
#' residuals, where residuals come from predicting outcomes and treatments with
#' covariates.
#'
#' @param data
#' @param y_vars
#' @param d_vars
#' @param x_vars
#' @param family
#' @param second_stage
#'
#' @return
#' @export
#'
#' @examples
residualise <- function(data, x_vars, y_var, family = 'ols'){
  x_data <- dplyr::select(data, {{x_vars}})
  y_data <- dplyr::pull(data, {{y_var}})

  df <- x_data %>%
    dplyr::mutate(y = y_data)

  # These following lines are stolen from Victor C.
  nobs <- nrow(x_data)
  foldid <- rep.int(1:2, times = ceiling(nobs / 2))[sample.int(nobs)]
  I <- split(1:nobs, foldid)
  fold1_indices <- I[[1]]
  fold2_indices <- I[[2]]


  y_hat <- ypreds <- rep(NA, nobs)

  thetas <- list()

  fold_1 <- df[fold1_indices, ]
  fold_2 <- df[fold2_indices, ] # Held-out Xs to predict values

  model_spec_fold1 <- NULL
  model_spec_fold2 <- NULL

  if(family == 'ols'){
    model_spec_fold1 <- parsnip::linear_reg() %>%
      parsnip::set_engine("lm")
    model_spec_fold2 <- model_spec_fold1

  } else if(family == 'rf'){

    model_spec_fold1 <- fit_cv_rf(fold_1)
    model_spec_fold2 <- fit_cv_rf(fold_2)

  }

  model_fit_fold1 <- parsnip::fit(model_spec_fold1, y ~ ., data = fold_1)
  model_fit_fold2 <- parsnip::fit(model_spec_fold2, y ~ ., data = fold_2)

  y_hat[fold2_indices] <- predict(model_fit_fold1, fold_2)$.pred - fold_2$y
  y_hat[fold1_indices] <- predict(model_fit_fold2, fold_1)$.pred - fold_1$y

  out <- list()

  out$model_fit_fold1 <- model_fit_fold1
  out$model_fit_fold2 <- model_fit_fold2
  out$y_hat <- y_hat
}

fit_cv_rf <- function(data
                      , grid_size = 25
                      , fold_size = 4){

  model_spec_tune <- parsnip::rand_forest(trees = tune(), min_n = tune()) %>%
    parsnip::set_mode("regression") %>%
    parsnip::set_engine("ranger")


  rf_cv <- rsample::vfold_cv(data, v = fold_size)
  #rf_grid <- expand.grid(trees = seq(500, 1000, length.out = 25), min_n = 5:10)

  dml_recipe <- recipes::recipe(y ~ ., data = data)

  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- parallel::makePSOCKcluster(all_cores)
  doParallel::registerDoParallel(cl)

  rf_res <- tune::tune_grid(object = model_spec_tune
                            , preprocessor = dml_recipe
                            , resamples = rf_cv
                            , grid = grid_size)

  best_hyper_parameters <- select_best(rf_res, metric = "rmse")

  model_spec <- parsnip::rand_forest(trees = best_hyper_parameters$trees, min_n = best_hyper_parameters$min_n) %>%
    parsnip::set_mode("regression") %>%
    parsnip::set_engine("ranger")

  model_spec
}
