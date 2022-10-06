
#
#' Estimates ATE and SE by regressing outcome residuals on treatment residuals
#'
#' @param y_resid
#' @param d_resid
#'
#' @return
#' @export
#'
#' @examples
dml_second_stage <- function(y_resid,
                             d_resid) {
  final_lm <- stats::lm(y_resid ~ d_resid)
  coef.est <- coef(final_lm)[2] # Extract coefficient
  se <- sqrt(sandwich::vcovHC(final_lm)[2, 2]) # Record standard error
  cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est, se))
  out <- list(
    coefficient = as.vector(coef.est),
    se = as.vector(se)
  )
  return(out)
}

#
#' This function facilitates running second-stage DML with LASSO
#'
#' @param formula
#' @param df
#' @param nfolds
#' @param lambda
#' @param alpha
#' @param standardize
#' @param post.lasso
#'
#' @return
#' @export
#'
#' @examples
dml_second_stage_lasso <- function(formula,
                                   df,
                                   nfolds,
                                   lambda = "1se",
                                   alpha = .05,
                                   standardize = TRUE,
                                   post.lasso = FALSE) {
  # Get sparse matrix
  x <- Matrix::sparse.model.matrix(formula, data = df)
  x <- x[, unique(Matrix::summary(x)$j)][, -1]
  # Print design matrix info
  cat(
    "\rDesign matrix has",
    dim(x)[1],
    "rows and",
    dim(x)[2],
    "columns\n"
  )
  # Run lasso
  set.seed(202103)
  final_lasso <- gamlr::cv.gamlr(
    x = x,
    y = df %>%
      dplyr::pull(obsDML:::formula_lhs(formula)),
    nfold = nfolds,
    standardize = standardize
  )
  # Get non-zero coefficienst
  c <- SparseM::as.matrix(coef(final_lasso, select = lambda))
  selected <- stats::setNames(
    c %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tibble::as_tibble(),
    c("variables", "coef")
  ) %>%
    dplyr::filter(coef != 0)
  if (post.lasso) {
    x_lm <- x[, selected$variables[-1]] %>%
      SparseM::as.matrix() %>%
      tibble::as_tibble()
    ols <- tryCatch(
      estimatr::lm_robust(df %>% dplyr::pull(obsDML:::formula_lhs(formula)) ~ ., data = x_lm),
      error = function(e) {
        cat("\rPost-LASSO estimator threw an error :(\t\t\t\n")
        return(NULL)
      }
    )
    if (is.null(ols)) {
      return(list(
        "lasso_selected" = selected,
        "lasso_model" = final_lasso
      ))
    }
    ols_out <- broom::tidy(ols, conf.int = TRUE, conf.level = 1 - alpha)
  }
  if (post.lasso) {
    return(list(
      "lasso_selected" = selected,
      "lasso_model" = final_lasso,
      "ols_model" = ols,
      "ols_coefficients" = ols_out
    ))
  } else {
    return(list(
      "lasso_selected" = selected,
      "lasso_model" = final_lasso
    ))
  }
}

#
#' This function facilitates running second-stage DML with Ridge
#'
#' @param formula
#' @param df
#' @param nfolds
#' @param lambda
#' @param alpha
#' @param standardize
#'
#' @return
#' @export
#'
#' @examples
dml_second_stage_ridge <- function(formula,
                                   df,
                                   nfolds,
                                   lambda = "lambda.1se",
                                   alpha = .05,
                                   standardize = TRUE) {
  # Get sparse matrix
  x <- Matrix::sparse.model.matrix(formula, data = df)
  x <- x[, unique(Matrix::summary(x)$j)][, -1]
  # Design matrix info
  cat(
    "\rDesign matrix has",
    dim(x)[1],
    "rows and",
    dim(x)[2],
    "columns\n"
  )
  # Run ridge
  set.seed(202103)
  final_ridge <- glmnet::cv.glmnet(
    x = x,
    y = df %>%
      dplyr::pull(obsDML:::formula_lhs(formula)),
    nfolds = nfolds,
    alpha = 0,
    standardize = standardize
  )
  # Get coefficients
  c <- SparseM::as.matrix(coef(final_ridge, s = lambda))
  selected <- stats::setNames(
    c %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tibble::as_tibble(),
    c("variables", "coef")
  )
  return(list(
    "ridge_selected" = selected,
    "ridge_model" = final_ridge
  ))
}


#' Runs bootstrapped hte_ridge to get coefficients and conf. intervals for them.
#'
#' @param formula
#' @param df
#' @param nfolds
#' @param lambda
#' @param alpha
#' @param nboot
#' @param standardize
#' @param strata
#'
#' @return
#' @export
#'
#' @examples
dml_second_stage_ridge_boot <- function(formula,
                                        df,
                                        nfolds,
                                        lambda = "lambda.min",
                                        alpha = .05,
                                        nboot = 1000,
                                        standardize = TRUE,
                                        strata) {
  # Get sparse matrix
  x <- Matrix::sparse.model.matrix(formula, data = df)
  x <- x[, unique(Matrix::summary(x)$j)][, -1]
  # Design matrix info
  cat(
    "\rDesign matrix has",
    dim(x)[1],
    "rows and",
    dim(x)[2],
    "columns\n"
  )
  y <- df %>% dplyr::pull(obsDML:::formula_lhs(formula))
  # Define function that runs ridge on bootstrap sample
  set.seed(202103)
  glmnet_call <- function(idx) {
    # idx <- sort(sample(1:nrow(x), replace = TRUE))
    boot_mod <- glmnet::cv.glmnet(
      x = x[idx, ],
      y = y[idx],
      alpha = 0,
      nfolds = 5,
      standardize = standardize
    )
    c <- SparseM::as.matrix(coef(boot_mod, s = lambda))
    selected <- stats::setNames(
      c %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tibble::as_tibble(),
      c("variables", "coef")
    )
    return(selected)
  }
  # Get Coefficients
  initial_model <- glmnet_call(1:nrow(x))
  # Run ridge regression across each bootstrap sample
  boot_models <- future.apply::future_lapply(1:nboot, function(i) {
    idx <- timeDate::sample(1:nrow(x), nrow(x), replace = TRUE)
    glmnet_call(idx)
  }, future.seed = TRUE)
  # Merge all n dataframes of ridge estimates
  # Every row is n estimates of a variable's coefficient
  boot_results <- boot_models %>%
    purrr::reduce(left_join, by = "variables") %>%
    dplyr::mutate(dplyr::across(-variables,
                  ~ .x - initial_model$coef))
  # Convert results into a matrix and calculate row quantiles
  # based on the user-specified value of alpha. Creates UB and LB
  LB <- alpha / 2
  UB <- 1 - LB
  row_quantiles <- boot_results %>%
    dplyr::select(-variables) %>%
    SparseM::as.matrix() %>%
    apply(1,
          quantile,
          c(1 - LB, 1 - UB),
          na.rm = TRUE
    ) %>%
    SparseM::t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(),
                  ~ initial_model$coef - .x)) %>%
    `names<-`(paste0(100*c(LB, UB), "%"))
  # Combine Bootstrapped confidence intervals with the estimated coefficients
  full_results <- dplyr::bind_cols(
    row_quantiles,
    coef = initial_model$coef,
    "variables" = boot_results$variables
  ) %>%
    `[`(, c(4, 1, 3, 2))
  return(full_results)
}





# Plot histogram of propensity scores
prop_scores_hist <- function(data){
  data %>%
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
}
# Plot Actual y values vs Predicted y values
actual_vs_pred <- function(data){
  data %>%
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
}
# Plot histogram of Outcome residuals
outcome_resid_hist <- function(data){
  data %>%
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
}
# Plot histogram of Treatment residuals
treat_resid_hist <- function(data){
  data %>%
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
}


