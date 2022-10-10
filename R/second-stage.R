



# dml(data, y_var, x_vars, d_vars, family = c('ols', ‘rf', ‘lasso’), second_stage = c(’mr', ‘sr1', 'sr2’))
#
# ‘mr’ - multiple residualsations y_hat ~ d_hat_i
#
# 'sr1' - single residualisation (old) y_hat ~d_hat_total:d_i
#
# ‘sr2’ - single residualisation (new) y_hat ~ d_i + d_hat_total
#
# returns
# first stage predictions
# lm object from stage 2
# tidy::broom version of stage 2 estimates

#
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
dml.lm <- function(data, y_var, x_vars, d_vars,
                   family = c('ols', 'rf', 'lasso'),
                   second_stage = c('mr', 'sr1', 'sr2')){
  if(second_stage == 'mr'){
    # Multiple Residualization
    # y_hat ~ sum_j d_hat_j
    yd_hat <- predict.dml(data, c(y_var, d_vars), x_vars, family=family)$y_hat
    reg.data = data[, c(y_var, d_vars)] - yd_hat
    colnames(reg.data) = c(y_var, d_vars)
    reg.str = paste('y ~ ', paste(d_vars, collapse=' + '), sep='')
  }else if(second.stage %in% c('sr1', 'sr2')){
    data.0 = data
    data.0$d_bar = apply(d_vars, MARGIN=1, FUN=sum)
    if(!all(d_bar %in% c(0, 1))){
      stop('dml.lm() 47: d_vars not partition')
    }
    yd_hat <- predict.dml(data.0, c(y_var, 'd_bar'), x_vars, family=family)$y_hat
    if(second_stage == 'sr1'){
      # Single Residualization - Old
      # y_hat ~ sum_j d_j*d_bar_hat
      reg.data = cbind(data.0[, y_var] - yd_hat[, 1:length(y_var)],
                       data.0$d_bar - yd_hat[, ncol(yd_hat)],
                       data.0[, d_vars] * replicate(length(d_vars), data.0$d_bar - yd_hat[, ncol(yd_hat)]))
      colnames(reg.data) = c(y_var, 'd_bar', d_vars)
      reg.str = paste('y ~ ', paste(d_vars, collapse=' + '), sep='')
    }else if(second_stage == 'sr2'){
      # Single Residualization - New
      # y_hat ~ d_bar_hat + sum_j d_j
      reg.data = cbind(data.0[, y_var] - yd_hat[, 1:length(y_var)],
                       data.0$d_bar - yd_hat[, ncol(yd_hat)],
                       data.0[, d_vars])
      colnames(reg.data) = c(y_var, 'd_bar', d_vars)
      reg.str = paste('y ~ d_bar + ', paste(d_vars, collapse=' + '), sep='')
    }else{
      stop('dml.lm() 62: sr1 or sr2 but neither sr1 nor sr2? what?')
    }
  }else{
    stop('dml.lm() 65: second_stage should be mr, sr1, or sr2')
  }

  out = list()
  for(i in 1:length(y_var)){
    reg.data$y = reg.data[, y_var[i]]
    out[[y_var[i]]] = lm(reg.str, data=reg.data)
  }
  # final_lm <- stats::lm(y_resid ~ d_resid)
  # coef.est <- coef(final_lm)[2] # Extract coefficient
  # se <- sqrt(sandwich::vcovHC(final_lm)[2, 2]) # Record standard error
  # cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est, se))
  # out <- list(
  #   coefficient = as.vector(coef.est),
  #   se = as.vector(se)
  # )
  return(out)
}


