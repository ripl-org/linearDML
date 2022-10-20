
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
                   , second_stage_family = c('mr', 'sr1', 'sr2')){
  out = list()
  if(second_stage_family == 'mr'){
    # Multiple Residualization
    # y_hat ~ sum_j d_hat_j

    # Predict outcomes and treatments
    y_resids <- predict.dml(data, y_var = y_var, x_vars = x_vars, family = first_stage_family)$resids

    d_resids <- lapply(d_vars, function(d_var) predict.dml(data, y_var = d_var, x_vars = x_vars, family = first_stage_family)$resids)
    names(d_resids) <- d_vars

    # Residualize outcomes and treatments
    reg.data <- data.frame(y_resids, d_resids)

    # Covariates include all (residualized) treatments

    #reg.str.0 = paste(d_vars, collapse=' + ')
  }else if(second_stage_family %in% c('sr1', 'sr2')){

    data$d_bar <- apply(data %>% dplyr::select(d_vars), MARGIN=1, FUN=sum)
    if(!all(data$d_bar %in% c(0, 1))){
      stop('dml.lm() 41: d_vars not partition')
    }

    # Predict outcomes and overall treatment (sum of treatments)
    #yd_hat <- residualise(data.0, c(y_var, 'd_bar'), x_vars, family=family)$y_hat
    #out$predicted.values = yd_hat

    y_resids <- predict.dml(data, y_var = y_var, x_vars = x_vars, family = first_stage_family)$resids
    d_bar_output <- predict.dml(data, y_var = 'd_bar', x_vars = x_vars, family = first_stage_family)
    d_bar_resids <- d_bar_output$resids
    d_bar_hat <- d_bar_output$predictions

    if(second_stage_family == 'sr1'){
      # Single Residualization - Old
      # y_hat ~ sum_j d_j*d_bar_hat

      # Residualize outcomes, and multiply treatments by residualized overall treatment
      #reg.data = cbind(data.0[, y_var] - yd_hat[, 1:length(y_var)],
      #                 data.0$d_bar - yd_hat[, ncol(yd_hat)],
      #                 data.0[, d_vars] * replicate(length(d_vars), data.0$d_bar - yd_hat[, ncol(yd_hat)]))

      reg.data <- data %>%
        dplyr::select(d_vars) %>%
        dplyr::mutate_all(~.x*d_bar_resids) %>% #interact all d's with d_bar_resids
        dplyr::mutate(y_resids = y_resids)


      # Covariates include all treatments (multiplied by overall residuals)
      #reg.str.0 = paste(d_vars, collapse=' + ')
    }else if(second_stage_family == 'sr2'){

      # Single Residualization - New
      # y_hat ~ d_bar_hat + sum_j d_j

      # Residualize outcomes
      #reg.data = cbind(data.0[, y_var] - yd_hat[, 1:length(y_var)],
      #                 yd_hat[, ncol(yd_hat)],
      #                 data.0[, d_vars])
      #colnames(reg.data) = c(y_var, 'd_bar', d_vars)

      reg.data <- data %>%
        dplyr::select(d_vars) %>%
        dplyr::mutate(d_bar_hat = d_bar_hat
               , y_resids = y_resids)

      # Covariates include all treatments and predicted overall treatment
      #reg.str.0 = paste('d_bar + ', paste(d_vars, collapse=' + '), sep='')
    }else{
      stop('dml.lm() 73: sr1 or sr2 but neither sr1 nor sr2? what?')
    }
  }else{
    stop('dml.lm() 76: second_stage should be mr, sr1, or sr2')
  }

  #out$model = list()
  #for(i in 1:length(y_var)){
    # Generate formula with current outcome and appropriately list of covariates
  #  reg.str.1 = paste(y_var[i], ' ~ ', reg.str.0, sep='')

    # Fit a linear model where formula is determined above by 'second_stage' type
   # out$model[[y_var[i]]] = lm(reg.str, data=reg.data)
    # Since we return the entire lm() object, no need to calculate specific
    # standard errors.
  #}

  model <- lm(data = reg.data, y_resids ~ .)
  return(model)
}



