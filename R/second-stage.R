
#' Estimates treatment effects by regressing outcome residuals on treatment
#' residuals, where residuals come from predicting outcomes and treatments with
#' covariates.
#'
#' @param data
#' @param y_vars
#' @param d_vars
#' @param x_vars
#' @param first_stage_family
#' @param second_stage_family
#' @param use_default_hyper
#'
#' @return
#' @export
#'
#' @examples
dml.lm <- function(data
                   , y_var
                   , d_vars
                   , x_vars
                   , first_stage_family = c('ols', 'rf', 'lasso')
                   , second_stage_family = c('mr', 'sr1', 'sr2')
                   , ...){
  out = list()

  y_model <- predict.dml(data, y_var = y_var, x_vars = x_vars, family = first_stage_family, ...)
  y_resids <- y_model$resids

  if(second_stage_family == 'mr'){
    # Multiple Residualization
    # y_hat ~ sum_j d_hat_j

    # Predict outcomes and treatments
    d_models <- lapply(d_vars, function(d_var) predict.dml(data
                                                           , y_var = d_var
                                                           , x_vars = x_vars
                                                           , family = first_stage_family
                                                           , ...))

    d_resids <- d_models %>% lapply(purrr::pluck, 'resids')
    names(d_resids) <- d_vars
    out$d_model <- d_models
    # Residualize outcomes and treatments
    reg.data <- data.frame(y_resids, d_resids)

    # Covariates include all (residualized) treatments

  }else if(second_stage_family %in% c('sr1', 'sr2')){

    data$d_bar <- apply(data %>% dplyr::select(d_vars), MARGIN=1, FUN=sum)
    if(!all(data$d_bar %in% c(0, 1))){
      stop('dml.lm(): d_vars not partition')
    }

    # Predict outcomes and overall treatment (sum of treatments)
    d_bar_output <- predict.dml(data
                                , y_var = 'd_bar'
                                , x_vars = x_vars
                                , family = first_stage_family
                                , ...)
    d_bar_resids <- d_bar_output$resids
    d_bar_hat <- d_bar_output$predictions
    out$d_model <- d_bar_output
    if(second_stage_family == 'sr1'){
      # Single Residualization - Old
      # y_hat ~ sum_j d_j*d_bar_hat

      # Residualize outcomes, and multiply treatments by residualized overall treatment
      reg.data <- data %>%
        dplyr::select(d_vars) %>%
        dplyr::mutate_all(~.x*d_bar_resids) %>% #interact all d's with d_bar_resids
        dplyr::mutate(y_resids = y_resids)


    }else if(second_stage_family == 'sr2'){

      # Single Residualization - New
      # y_hat ~ d_bar_hat + sum_j d_j
      # Covariates include all treatments and predicted overall treatment
      reg.data <- data %>%
        dplyr::select(d_vars) %>%
        dplyr::mutate(d_bar_hat = d_bar_hat
                      , y_resids = y_resids)


    }else{
      stop('dml.lm(): sr1 or sr2 but neither sr1 nor sr2? what?')
    }
  }else{
    stop('dml.lm(): second_stage should be mr, sr1, or sr2')
  }

  model <- lm(data = reg.data, y_resids ~ .)
  out$y_model <- y_model
  out$model <- model

  return(out)
}



