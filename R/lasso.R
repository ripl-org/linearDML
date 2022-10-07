#
#
#' Random forest w/ grid-search possibility
#' Grid search ranger
#' @param form
#' @param df
#' @param probability
#' @param predict_df
#' @param mtry
#' @param node_size


#' Title
#'
#' @param x
#' @param y
#' @param d
#' @param y_family
#' @param d_family
#' @param nfold
#' @param formula
#' @param free_vars
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
dml_first_stage_lasso <- function(data,
                                  x_vars,
                                  y_var,
                                  d_vars,
                                  y_family,
                                  d_family,
                                  nfold = 2,
                                  verbose = TRUE,
                                  ...) {
  # Create model matrix. Removes constant columns and duplicates.
  x_data <- data %>% dplyr::select({{x_vars}})
  x_data <- Matrix::sparse.model.matrix(~ ., data = x_data)
  x_data <- x_data[, unique(Matrix::summary(x_data)$j)][, -1]
  y_data <- data %>% dplyr::pull({{y_var}})
  d_data <- data %>% dplyr::select({{d_vars}})


  if (verbose) {
    cat(
      "\rDesign matrix has",
      dim(x_data)[1],
      "rows and",
      dim(x_data)[2],
      "columns\n"
    )
  }
  # These following lines are stolen from Victor C.
  nobs <- nrow(x_data)
  foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)]
  I <- split(1:nobs, foldid)
  ytil <- ypreds <- rep(NA, nobs)
  dtil <- dpreds <- matrix(nrow = nobs, ncol = ncol(d_data))
  colnames(dtil) <- colnames(d_data)
  colnames(dpreds) <- colnames(d_data)
  thetas <- list()


  # Initialize dataframes for returning values AND selected predictors
  y_predictors <- d_predictors <- tibble::tibble()
  # The following loop builds residuals over the folds
  for (b in 1:length(I)) {

    # Following four lines subsets each component of data
    x_sub <- x_data[-I[[b]], ]
    y_sub <- y_data[-I[[b]]]
    d_sub <- d_data[-I[[b]], ]

    x_test <- x_data[I[[b]], ] # Held-out Xs to predict values
    y_test <- y_data[I[[b]]]
    d_test <- d_data[I[[b]], ]

    # Residualize
    y_model <- obsDML:::lasso_helper(
      x = x_sub,
      y = y_sub,
      family = y_family,
      test_x = x_test,
      ...
    )


    d_models <- map(1:ncol(d_sub),function(d_index){

      d_model <- obsDML:::lasso_helper(
        x = x_sub,
        y = d_sub[,d_index],
        family = d_family,
        test_x = x_test,
        ...
      )

    })

    # Return predicted values
    dhat <- d_models %>% sapply(pluck, 'values')
    yhat <- y_model$values

    # Return selected LASSO variables
    d_predictors <- dplyr::bind_rows(d_predictors, d_models %>% map_df(pluck, 'coefficients'))
    y_predictors <- dplyr::bind_rows(y_predictors, y_model$coefficients)

     # Get predicted values
    ypreds[I[[b]]] <- yhat
    dpreds[I[[b]],] <- dhat

    # Calculate Residuals
    ytil[I[[b]]] <- (y_test - yhat)
    dtil[I[[b]],] <- (as.matrix(d_test) - dhat)

    partial_model_df <- tibble::tibble(y_til = ytil[I[[b]]]) %>%
      dplyr::bind_cols(dtil[I[[b]],])

    partial_model_coefs <- stats::lm(data = partial_model_df,  y_til ~ .) %>%
      summary %>%
      coef
    # Calculate theta for DML1 calculation
    theta <- list(
      theta = partial_model_coefs[, 1][-1],
      se = partial_model_coefs[, 2][-1],
      ytil = ytil[I[[b]]],
      dtil = dtil[I[[b]],]
    )
    thetas <- append(thetas, list("theta" = theta))
    if (verbose == TRUE) {
      cat("\r--- Finished fold", b, "---\n")
    }
  }

  browser()
  # List of error-checking
  obsDML:::assert(sum(is.na(ytil)) == 0)
  obsDML:::assert(sum(is.na(dtil)) == 0)
  obsDML:::assert(sum(is.na(ypreds)) == 0)
  obsDML:::assert(sum(is.na(dpreds)) == 0)
  obsDML:::assert(sum(is.na(thetas)) == 0, message = "Vector 'thetas' has NA values")
  obsDML:::assert(length(ytil) == nobs)
  # Calculate precision metrics

  outcomes <- d_data %>%
    dplyr::rename_all(~ paste(., "truth", sep = "_")) %>%
    dplyr::bind_cols(tibble::as_tibble(dpreds) %>% dplyr::rename_all(~ paste(., "pred", sep = "_"))) %>%
    dplyr::bind_cols(tibble::as_tibble(dtil) %>% dplyr::rename_all(~ paste(., "resid", sep = "_"))) %>%
    dplyr::mutate(observation = dplyr::row_number()
                  , y_truth = y_data
                  , y_resid = ytil
                  , y_pred = ypreds)


  rsquared <- yardstick::rsq(
    data = outcomes,
    y_truth,
    y_pred
  )$`.estimate`

  RMSE <- yardstick::rmse(
    data = outcomes,
    y_truth,
    y_pred
  )$`.estimate`

  #outcome_resid_hist <- outcome_resid_hist(y_outcomes)

  d_outcomes <- outcomes %>%
    tidyr::pivot_longer(c(-observation, -y_truth, -y_resid, -y_pred), names_sep = '_', names_to = c('name', 'type')) %>%
    dplyr::group_by(name) %>%
    tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data)
      data %>%
        tidyr::pivot_wider(id_cols = c(observation, y_truth, y_resid, y_pred), names_from = type, values_from = value) %>%
        dplyr::mutate(d_truth = truth
                      , d_pred = pred))) %>%
    dplyr::mutate(auc = purrr::map_dbl(data, function(data)obsDML:::auc_roc(data$d_pred, data$d_truth))
                  , outcome_resid_hist = purrr::map(data,outcome_resid_hist)
                  , prop_scores_hist = purrr::map(data, prop_scores_hist)
                  , actual_vs_pred = purrr::map(data, actual_vs_pred)
                  , treat_resid_hist = purrr::map(data, treat_resid_hist))



  # Return residuals and predicted values
  out <- list(
    y_predicted_values = ypreds,
    d_predicted_values = dpreds,
    y_residuals = ytil,
    d_residuals = dtil,
    y_actual = outcomes$y_truth,
    d_actual = d_data,
    R2 = rsquared,
    RMSE = RMSE,
    d_propensity_hist = prop_scores_hist,
    d_outcomes = d_outcomes,
    theta_estimates = thetas
  )
}




#' Define function that runs k-fold validated lasso
#'
#' @param x
#' @param y
#' @param nfolds
#' @param family
#' @param test_x
#' @param free_idx
#' @param stand
#'
#' @return
#'
#' @examples
lasso_helper <- function(x,
                         y,
                         nfolds = 5,
                         family,
                         test_x = NULL,
                         ...) {
  set.seed(202103)
  model <- gamlr::cv.gamlr(
    x = x,
    y = y,
    nfold = nfolds,
    family = family,
  ) # Specifies variables to keep
  # Predict values using lambda value that minimizes cv error
  if (!is.null(test_x)) {
    if (family == "binomial") {
      values <- predict(model, test_x, select = "1se", type = "response")
    } else {
      values <- predict(model, test_x, select = "1se")
    }
  }
  # Uncomment below to return both values and selected variables
  # Generate tibble for non-zero coefficients
  non_zero_coef <- tibble::tibble()
  # Get non-zero coefficienst
  c <- SparseM::as.matrix(coef(model, select = "1se"))
  non_zero_coef <- dplyr::bind_rows(
    non_zero_coef,
    tryCatch(
      tibble::tibble(
        selected_variables =
          names(c[Matrix::which(c != 0), ]),
        selected_coefs =
          c[Matrix::which(c != 0), ]
      ),
      error = function(e) {
        tibble::tibble(selected_variables = "None")
      }
    )
  )
  # Return predicted values
  out <- list(
    coefficients = non_zero_coef,
    model = model
  )
  if (!is.null(test_x)) {
    out <- append(out, list(values = as.vector(values)))
  }
  return(out)
}
