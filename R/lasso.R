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
  ytil <- ypreds <- dtil <- dpreds <- rep(NA, nobs)
  thetas <- list()


  # Initialize dataframes for returning values AND selected predictors
  y_predictors <- d_predictors <- tibble::tibble()
  # The following loop builds residuals over the folds
  for (b in 1:length(I)) {
    browser()
    # Following four lines subsets each component of data
    x_sub <- x_data[-I[[b]], ]
    y_sub <- y_data[-I[[b]]]
    d_sub <- d_data[-I[[b]], ]
    test_x <- x_data[I[[b]], ] # Held-out Xs to predict values
    # Residualize
    y_model <- obsDML:::lasso_helper(
      x = x_sub,
      y = y_sub,
      family = y_family,
      test_x = test_x,
      ...
    )



    dhat <- sapply(1:ncol(test_d),function(d_index){

      d_model <- obsDML:::lasso_helper(
        x = x_sub,
        y = d_sub,
        family = d_family,
        test_x = test_x,
        ...
      )


      d_model$values
    })

    # Return predicted values
    yhat <- y_model$values
    # Return selected LASSO variables
    d_predictors <- dplyr::bind_rows(d_predictors, d_model$coefficients)
    y_predictors <- dplyr::bind_rows(y_predictors, y_model$coefficients)
    # Get predicted values
    ypreds[I[[b]]] <- yhat
    dpreds[I[[b]]] <- dhat
    # Calculate Residuals
    ytil[I[[b]]] <- (obsDML:::as_numeric(y[I[[b]]]) - yhat)
    dtil[I[[b]]] <- (obsDML:::as_numeric(d[I[[b]]]) - dhat)
    # Calculate theta for DML1 calculation
    theta <- list(
      theta = coef(Matrix::summary(stats::lm(ytil[I[[b]]] ~ dtil[I[[b]]])))[, 1][2],
      se = coef(Matrix::summary(stats::lm(ytil[I[[b]]] ~ dtil[I[[b]]])))[, 2][2],
      ytil = ytil[I[[b]]],
      dtil = dtil[I[[b]]]
    )
    thetas <- append(thetas, list("theta" = theta))
    if (verbose == TRUE) {
      cat("\r--- Finished fold", b, "---\n")
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
  outcomes <- tibble::tibble(
    y_truth = y,
    y_pred = ypreds,
    d_truth = d,
    `0` = 1 - dpreds,
    d_pred = dpreds,
    y_resid = ytil,
    d_resid = dtil
  )
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
  auc <- obsDML::auc_roc(preds = outcomes$d_pred,
                         actuals = obsDML:::as_numeric(outcomes$d_truth))
  # Plot histogram of propensity scores
  prop_scores_hist <- outcomes %>%
    dplyr::rename("Treatment" = d_truth) %>%
    ggplot2::ggplot(ggplot2::aes(x = d_pred, fill = Treatment)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::facet_wrap(~Treatment, nrow = 2) +
    ggplot2::labs(
      x = "Propensity Scores",
      y = "Frequency",
      title = "Histogram of Propensity Scores"
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "gray95"),
      legend.key = ggplot2::element_rect(fill = "white"),
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  # Plot Actual y values vs Predicted y values
  actual_vs_pred <- outcomes %>%
    dplyr::rename("Treatment" = d_truth) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = y_truth,
      y = y_pred,
      color = Treatment
    )) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::facet_wrap(~Treatment, nrow = 2) +
    ggplot2::labs(
      x = "Y Truth",
      y = "Y Predicted",
      title = "Actual Y values vs Predicted Y values"
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "gray95"),
      legend.key = ggplot2::element_rect(fill = "white"),
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  # Plot histogram of Outcome residuals
  outcome_resid_hist <- outcomes %>%
    dplyr::rename("Treatment" = d_truth) %>%
    ggplot2::ggplot(ggplot2::aes(x = y_resid, fill = Treatment)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::facet_wrap(~Treatment, nrow = 2) +
    ggplot2::labs(
      x = "Outcome Residuals",
      y = "Frequency",
      title = "Histogram of Outcome Residuals"
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "gray95"),
      legend.key = ggplot2::element_rect(fill = "white"),
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  # Plot histogram of Treatment residuals
  treat_resid_hist <- outcomes %>%
    dplyr::rename("Treatment" = d_truth) %>%
    ggplot2::ggplot(ggplot2::aes(x = d_resid, fill = Treatment)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::facet_wrap(~Treatment, nrow = 2) +
    ggplot2::labs(
      x = "Treatment Residuals",
      y = "Frequency",
      title = "Histogram of Treatment Residuals"
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "gray95"),
      legend.key = ggplot2::element_rect(fill = "white"),
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  # Return residuals and predicted values
  out <- list(
    y_predicted_values = ypreds,
    d_predicted_values = dpreds,
    y_residuals = ytil,
    d_residuals = dtil,
    y_actual = outcomes$y_truth,
    d_actual = outcomes$d_truth,
    R2 = rsquared,
    RMSE = RMSE,
    AUC = auc,
    d_propensity_hist = prop_scores_hist,
    y_actual_vs_pred = actual_vs_pred,
    y_resid_hist = outcome_resid_hist,
    d_resid_hist = treat_resid_hist,
    y_selected_vars = suppressMessages(
      y_predictors %>%
        dplyr::group_by(selected_variables) %>%
        dplyr::summarise(
          n = n(),
          coef = Matrix::mean(selected_coefs)
        ) %>%
        dplyr::ungroup()
    ),
    d_selected_vars = suppressMessages(
      d_predictors %>%
        dplyr::group_by(selected_variables) %>%
        dplyr::summarise(
          n = n(),
          coef = Matrix::mean(selected_coefs)
        ) %>%
        dplyr::ungroup()
    ),
    theta_estimates = thetas
  )
  return(out)
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
