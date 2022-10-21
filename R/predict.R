
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
#' @param use_default_hyper
#' @return
#' @export
#'
#' @examples
predict.dml <- function(data
                        , y_var
                        , x_vars
                        , family = 'ols'
                        , ...){
  x_data <- dplyr::select(data, x_vars)
  y_data <- dplyr::pull(data, y_var)

  if(!family %in% c('ols', 'rf'))
    stop('argument "family" must be one of "ols" or "rf')

  if(is.logical(y_data)){
    warning(paste0('Logical variable "', y_var, '" being converted to 1/0 numeric'))
    y_data <- as.numeric(y_data)
  }

  df <- x_data %>%
    dplyr::mutate(y = y_data)

  #split data into two folds
  nobs <- nrow(x_data)
  foldid <- rep.int(1:2, times = ceiling(nobs / 2))[sample.int(nobs)]
  I <- split(1:nobs, foldid)
  fold1_indices <- I[[1]]
  fold2_indices <- I[[2]]


  predictions <- rep(NA, nobs)
  resids <- rep(NA, nobs)

  fold_1 <- df[fold1_indices, ]
  fold_2 <- df[fold2_indices, ] # Held-out Xs to predict values

  model_spec_fold1 <- NULL
  model_spec_fold2 <- NULL

  #fit model family
  if(family == 'ols'){
    model_spec_fold1 <- parsnip::linear_reg() %>%
      parsnip::set_engine("lm")

    model_spec_fold2 <- model_spec_fold1

  } else if(family == 'rf'){

    #perform cross validation on rf
    model_spec_fold1 <- fit_cv_rf(fold_1, ...)
    model_spec_fold2 <- fit_cv_rf(fold_2, ...)

  }

  #use the parameters of fold1 on fold2 data
  model_fit_fold1 <- parsnip::fit(model_spec_fold1, y ~ ., data = fold_1)
  model_fit_fold2 <- parsnip::fit(model_spec_fold2, y ~ ., data = fold_2)

  predictions[fold2_indices] <- predict(model_fit_fold1, fold_2)$.pred
  predictions[fold1_indices] <- predict(model_fit_fold2, fold_1)$.pred

  resids[fold2_indices] <- predictions[fold2_indices] - y_data[fold2_indices]
  resids[fold1_indices] <- predictions[fold1_indices] - y_data[fold1_indices]

  #save model specs and residuals
  out <- list(model_fit_fold1 = model_fit_fold1
              , model_fit_fold2 = model_fit_fold2
              , predictions = predictions
              , resids = resids)

  out
}

fit_cv_rf <- function(data
                      , grid_size = 20
                      , fold_size = 4
                      , use_default_hyper = TRUE){

  model_spec_default <- parsnip::rand_forest() %>%
    parsnip::set_mode("regression") %>%
    parsnip::set_engine("ranger")

  if(use_default_hyper)
    return(model_spec_default)

  message("Tuning random forest hyperparameters")
  model_spec_tune <- parsnip::rand_forest(trees = tune(), min_n = tune()) %>%
    parsnip::set_mode("regression") %>%
    parsnip::set_engine("ranger")


  rf_cv <- rsample::vfold_cv(data, v = fold_size)

  dml_recipe <- recipes::recipe(y ~ ., data = data)

  #all_cores <- parallel::detectCores(logical = FALSE)
  #cl <- parallel::makePSOCKcluster(all_cores)
  #doParallel::registerDoParallel(cl)

  rf_res <- tune::tune_grid(object = model_spec_tune
                            , preprocessor = dml_recipe
                            , resamples = rf_cv
                            , grid = grid_size
                            , control = control_grid(verbose = TRUE))

  best_hyper_parameters <- select_best(rf_res, metric = "rmse")

  model_spec <- parsnip::rand_forest(trees = best_hyper_parameters$trees, min_n = best_hyper_parameters$min_n) %>%
    parsnip::set_mode("regression") %>%
    parsnip::set_engine("ranger")

  model_spec
}
