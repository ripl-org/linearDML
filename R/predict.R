
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
                        , seed = 10272022
                        , ...){
  x_data <- dplyr::select(data, dplyr::all_of(x_vars))
  y_data <- dplyr::pull(data, y_var)
  mode <- 'regression'
  prediction_type <- 'numeric'

  if(!family %in% c('ols', 'rf'))
    stop('argument "family" must be one of "ols" or "rf')

  if(is.logical(y_data)){
    warning(paste0('Logical variable "', y_var, '" being converted to factor'))
    y_data <- as.factor(y_data)
    mode <- 'classification'
    prediction_type <- 'prob'
  }

  if('character' %in% lapply(x_data, class)){
    warning(paste0('Converting character columns to factor'))
    x_data <- x_data %>%
      dplyr::mutate_if(is.character, as.factor)
  }

  set.seed(seed)

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
    if(mode == 'classification')
      model_spec_fold1 <- parsnip::logistic_reg() %>%
        parsnip::set_engine("glm") %>%
        parsnip::set_mode(mode)

    if(mode == 'regression')
      model_spec_fold1 <- parsnip::linear_reg() %>%
        parsnip::set_engine("lm") %>%
        parsnip::set_mode(mode)

    model_spec_fold2 <- model_spec_fold1

  } else if(family == 'rf'){

    #perform cross validation on rf
    model_spec_fold1 <- fit_cv_rf(fold_1, rf_mode = mode, ...)
    model_spec_fold2 <- fit_cv_rf(fold_2, rf_mode = mode, ...)

  }

  #use the parameters of fold1 on fold2 data
  model_fit_fold1 <- parsnip::fit(model_spec_fold1, y ~ ., data = fold_1)
  model_fit_fold2 <- parsnip::fit(model_spec_fold2, y ~ ., data = fold_2)
  if(mode == 'classification'){
    predictions[fold2_indices] <- predict(model_fit_fold1, fold_2, type = prediction_type)$.pred_TRUE
    predictions[fold1_indices] <- predict(model_fit_fold2, fold_1, type = prediction_type)$.pred_TRUE

    resids[fold2_indices] <- as.logical(y_data[fold2_indices]) - predictions[fold2_indices]
    resids[fold1_indices] <- as.logical(y_data[fold1_indices]) - predictions[fold1_indices]
  }else{

    predictions[fold2_indices] <- predict(model_fit_fold1, fold_2, type = prediction_type)$.pred
    predictions[fold1_indices] <- predict(model_fit_fold2, fold_1, type = prediction_type)$.pred

    resids[fold2_indices] <- y_data[fold2_indices] - predictions[fold2_indices]
    resids[fold1_indices] <- y_data[fold1_indices] - predictions[fold1_indices]
  }



  #save model specs and residuals
  out <- list(model_fit_fold1 = model_fit_fold1
              , model_fit_fold2 = model_fit_fold2
              , predictions = predictions
              , resids = resids
              , data = df)

  out
}

fit_cv_rf <- function(data
                      , rf_mode = "regression"
                      , use_default_hyper = TRUE
                      , trees = NULL
                      , min_n = NULL){


  model_spec_default <- parsnip::rand_forest() %>%
    parsnip::set_mode(mode = rf_mode) %>%
    parsnip::set_engine("ranger")

  if(use_default_hyper & is.null(trees) & is.null(min_n)){
    message("Using default random forest hyperparameters")
    return(model_spec_default)
  }

  if(!is.null(trees)&!is.null(min_n)){
    message("Using user supplied random forest hyperparameters")
    model_spec <- parsnip::rand_forest(trees = trees, min_n = min_n) %>%
      parsnip::set_mode(mode = rf_mode) %>%
      parsnip::set_engine("ranger")

    return(model_spec)
  }


  model_spec_tune <- parsnip::rand_forest(trees = tune(), min_n = tune()) %>%
    parsnip::set_mode(mode = rf_mode) %>%
    parsnip::set_engine("ranger")

  if(!is.null(trees))
    model_spec_tune  <- model_spec_tune %>% set_args(trees = trees)
  if(!is.null(min_n))
    model_spec_tune  <- model_spec_tune %>% update(min_n = min_n)


  message("Tuning random forest hyperparameters")


  best_hyper_parameters <- cross_validate(data, model_spec_tune, rf_mode)

  model_spec <- parsnip::rand_forest() %>%
    parsnip::set_mode(mode = rf_mode) %>%
    parsnip::set_engine("ranger")

  if(is.null(trees))
    model_spec <- model_spec %>% set_args(trees = best_hyper_parameters$trees)
  else
    model_spec <- model_spec %>% set_args(trees = trees)

  if(is.null(min_n))
    model_spec <- model_spec %>% set_args(min_n = best_hyper_parameters$min_n)
  else
    model_spec <- model_spec %>% set_args(min_n = min_n)

  model_spec
}

cross_validate <- function(data
                           , model_spec_tune
                           , mode
                           , grid_size = 4
                           , fold_size = 4){


  rf_cv <- rsample::vfold_cv(data, v = fold_size)

  dml_recipe <- recipes::recipe(y ~ ., data = data)

  #all_cores <- parallel::detectCores(logical = FALSE)
  #cl <- parallel::makePSOCKcluster(all_cores)
  #doParallel::registerDoParallel(cl)

  tuned_models <- tune::tune_grid(object = model_spec_tune
                                  , preprocessor = dml_recipe
                                  , resamples = rf_cv
                                  , grid = grid_size
                                  , control = tune::control_grid(verbose = TRUE))
  metric <- 'rmse'
  if(mode == 'regression')
    metric <- 'rmse'
  if(mode == 'classification')
    metric <- 'accuracy'

  best_hyper_parameters <- tune::select_best(tuned_models, metric = metric)
}
