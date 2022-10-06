#
#' Random forest first stage
#'
#' @param data
#' @param x
#' @param y
#' @param d
#' @param verbose
#' @param nfold
#' @param mtry
#' @param num_trees
#' @param node_size
#' @param error_type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
dml_first_stage_rf <- function(data,
                               x_vars,
                               y_var,
                               d_vars,
                               verbose = TRUE,
                               nfold = 2,
                               mtry = NULL,
                               num_trees= NULL,
                               node_size= NULL,
                               error_type = "OOB",
                               ...) {

  df <- dplyr::select(data, {{x_vars}} , y := {{y_var}})
  d_df <- dplyr::select(data, {{d_vars}})
  # These following lines are stolen from Victor C.
  nobs <- nrow(df)
  foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)]
  I <- split(1:nobs, foldid)
  ytil <- ypreds <- rep(NA, nobs)
  dtil <- dpreds <- matrix(nrow = nobs, ncol = ncol(d_df))
  colnames(dtil) <- colnames(d_df)
  colnames(dpreds) <- colnames(d_df)
  thetas <- list()
  yprob <- is.factor(df$y)

  # The following loop builds residuals over the folds
  for (b in 1:length(I)) {
    # Following four lines subsets each component of data
    x_sub <- df[-I[[b]], ]
    test_x <- df[I[[b]], ] # Held-out Xs to predict values

    d_sub <- d_df[-I[[b]], ]
    test_d <- d_df[I[[b]], ]

    # Residualize
    y_model <- suppressMessages(
      rf_helper(
        form = y ~ .,
        df = x_sub ,
        probability = yprob,
        mtry = mtry,
        num_trees = num_trees,
        node_size = node_size,
        predict_df = test_x ,
        error_type = error_type,
        verbose = FALSE
      )
    )

    yhat <- y_model$values
    dhat <- sapply(1:ncol(test_d),function(d_index){

      d_model_df <- x_sub %>%
        dplyr::select(-y) %>%
        dplyr::mutate(d = dplyr::pull(d_sub, d_index))

      d_predict_df <- test_x %>%
        dplyr::select(-y) %>%
        dplyr::mutate(d = dplyr::pull(test_d, d_index))

      dprob <- is.factor(d_model_df)

      d_model <- suppressMessages(
        rf_helper(
          form = d ~ .,
          df = d_model_df ,
          probability = dprob,
          mtry = mtry,
          num_trees = num_trees,
          node_size = node_size,
          predict_df = d_predict_df,
          error_type = error_type,
          verbose = FALSE
        )
      )

      d_model$values
    })

    # Get predicted values
    ypreds[I[[b]]] <- yhat
    dpreds[I[[b]],] <- dhat
    # Calculate Residuals
    ytil[I[[b]]] <- (obsDML:::as_numeric(test_x$y) - yhat)
    dtil[I[[b]],] <- (SparseM::as.matrix(test_d) - dhat)


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
      cat("--- Finished fold", b, "---\n")
    }
  }

  # List of error-checking
  obsDML:::assert(sum(is.na(ytil)) == 0)
  obsDML:::assert(sum(is.na(dtil)) == 0)
  obsDML:::assert(sum(is.na(ypreds)) == 0)
  obsDML:::assert(sum(is.na(dpreds)) == 0)
  obsDML:::assert(sum(is.na(thetas)) == 0, message = "Vector 'thetas' has NA values")
  obsDML:::assert(length(ytil) == nobs)
  # Calculate precision metrics

  outcomes <- d_df %>%
    dplyr::rename_all(~ paste(., "truth", sep = "_")) %>%
    dplyr::bind_cols(tibble::as_tibble(dpreds) %>% dplyr::rename_all(~ paste(., "pred", sep = "_"))) %>%
    dplyr::bind_cols(tibble::as_tibble(dtil) %>% dplyr::rename_all(~ paste(., "resid", sep = "_"))) %>%
    dplyr::mutate(observation = dplyr::row_number()
                  , y_truth = df$y
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
    d_actual = d_df,
    R2 = rsquared,
    RMSE = RMSE,
    d_propensity_hist = prop_scores_hist,
    d_outcomes = d_outcomes,
    theta_estimates = thetas
  )
  return(out)
}

#' @param num_trees
#' @param error_type
#' @param verbose
#'
#' @return
#'
#' @examples
rf_helper <- function(form,
                      df,
                      probability,
                      predict_df = NULL,
                      mtry = NULL,
                      node_size = NULL,
                      num_trees = NULL,
                      error_type = "OOB",
                      verbose = FALSE,
                      ...) {
  # Error checking
  obsDML:::assert(
    is.null(error_type) | error_type == "OOB",
    "Argument 'error_type' must be either NULL or 'OOB'"
  )
  # Split out data by train and test set
  x <- df %>% dplyr::select(-!!obsDML:::formula_lhs(form))
  y <- df %>% dplyr::pull(obsDML:::formula_lhs(form))
  if (!is.null(predict_df)) {
    predict_x <- predict_df %>% dplyr::select(-!!obsDML:::formula_lhs(form))
    predict_y <- predict_df %>% dplyr::pull(obsDML:::formula_lhs(form))
  }
  # Run models across tuning grid if specified
  # If not, run model with default values
  if (!is.null(mtry) | !is.null(node_size) | !is.null(num_trees)) {
    if (is.null(mtry)) {
      mtry <- floor(sqrt(ncol(x)))
    }
    if (is.null(node_size)) {
      node_size <- ifelse(probability == TRUE, 1, 5)
    }
    if (is.null(num_trees)) {
      num_trees <- 1000
    }
    # Define grid of hyper parameters
    grid <- SparseM::t(expand.grid(
      mtry = mtry,
      node_size = node_size,
      num_trees = num_trees
    )) %>%
      as.data.frame()
    # Run a random forest for every row in grid
    set.seed(202103)
    models <- future.apply::future_lapply(
      grid,
      function(i) {
        model <- ranger::ranger(y ~ .,
                                data = x,
                                probability = probability,
                                mtry = i[1],
                                min.node.size = i[2],
                                num.trees = i[3],
                                verbose = FALSE,
                                ...
        )
        error <- model$prediction.error
        return(error)
      }
    )
    # Append the respective OOB rate to each row of the grid
    grid_errors <- tibble::as_tibble(SparseM::t(grid)) %>%
      dplyr::bind_cols(error = unlist(models)) %>%
      dplyr::arrange(error)
    # Get the hyper parameters that generated minimum OOB error
    hyper_params <- tibble::as_tibble(SparseM::t(grid))[which.min(unlist(models)), ]
    # Run model with hyper params defined above
    final_model <- ranger::ranger(y ~ .,
                                  data = x,
                                  probability = probability,
                                  mtry = hyper_params$mtry,
                                  min.node.size = hyper_params$node_size,
                                  num.trees = hyper_params$num_trees
    )
  } else {
    if (verbose == TRUE) {
      cat(
        "At least one of mtry, node_size, and num_trees are null,",
        "so using default ranger values\n"
      )
    }
    # If hyper parameters aren't defined, run model with default values
    final_model <- ranger::ranger(y ~ .,
                                  data = x,
                                  probability = probability
    )
  }
  # If a dataframe is provided, generate predictions
  if (!is.null(predict_df)) {
    values <- if (is.factor(y)) {
      predict(final_model, predict_x)$predictions[, 2]
    } else {
      predict(final_model, predict_x)$predictions
    }
  }
  # List of outs
  out <- list(model = final_model)
  if (!is.null(predict_df)) out <- append(out, list(values = as.vector(values)))
  if (!is.null(mtry) & !is.null(node_size) & !is.null(num_trees)) {
    out <- append(
      out,
      list(
        grid = grid_errors,
        which_min = hyper_params
      )
    )
  }
  return(out)
}
