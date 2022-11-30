
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
                   , x_vars
                   , d_vars = NULL
                   , h_vars = NULL
                   , first_stage_family = c('ols', 'rf', 'lasso', 'user-defined')
                   , predict_fun = NULL
                   , second_stage_family = c('mr', 'sr1', 'sr2')
                   , ...){
  out = list()


  if(is.null(d_vars)){
    if(second_stage_family %in% c('sr1', 'sr2')){
      stop('dml.lm: if second_stage_family is in (sr1, sr2), d_vars cannot be NULL')
    }
    if(is.null(h_vars)){
      stop('dml.lm: if d_vars and h_vars cannot both be NULL')
    }
  }

  if(first_stage_family == 'user-defined'){
    if(is.null(predict_fun)){
      stop(paste('dml.lm: if prediction is user-defined, predict_fun must be',
                 'function(data, y_var, x_vars) returning list(predictions, resids)', sep=''))
    }
  }else{
    predict_fun <- function(data, y_var, x_vars){
      return(predict.dml(data, y_var = y_var, x_vars = x_vars, family = first_stage_family, ...))
    }
  }
  y_model <- predict_fun(data, y_var = y_var, x_vars = x_vars, ...)
  y_resids <- y_model$resids


  if(second_stage_family == 'mr'){
    if(is.null(h_vars)){
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

    }else{
      # Multiple Residualization, with covariate interactions
      # y_hat ~ sum_j h_j(x) d_hat_j

      d_vars = unique(h_vars$d)

      # Predict treatments
      d_models <- lapply(d_vars, function(d_var) predict.dml(data
                                                             , y_var = d_var
                                                             , x_vars = x_vars
                                                             , family = first_stage_family
                                                             , ...))

      d_resids <- d_models %>% lapply(purrr::pluck, 'resids')
      names(d_resids) <- d_vars
      out$d_model <- d_models

            #
      # h_vars should be a data frame with three columns:
      # d - the treatment variable
      # fx - the (transformation of the) covariates
      # fxd.name - the name of this combined variable in the second stage
      reg.data <- data.frame(y_resids)
      reg.data[, h_vars$fxd.name] <- data[, h_vars$fx] * purrr::pluck(d_resids, h_vars$d)
    }
  }else if(second_stage_family %in% c('sr1', 'sr2')){
    if(!is.null(h_vars)){
      stop('dml.lm(): treatment heterogeneity not defined for single residualization')
    }

    data$d_bar <- apply(data %>% dplyr::select(d_vars), MARGIN=1, FUN=sum)
    if(!all(data$d_bar %in% c(0, 1))){
      stop('dml.lm(): d_vars not partition')
    }

    # Predict outcomes and overall treatment (sum of treatments)
    d_bar_output <- predict_fun(data
                                , y_var = 'd_bar'
                                , x_vars = x_vars
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

  model <- lm(data = reg.data, y_resids ~ . - 1)
  out$y_model <- y_model
  out$model <- model

  return(out)
}


#' Estimates treatment effects by regressing outcome residuals on treatment
#' residuals, where residuals come from predicting outcomes and treatments with
#' covariates.
#'
#' @param formula
#' @param data
#' @param first_stage_family
#' @param second_stage_family
#' @param use_default_hyper
#'
#' @return
#' @export
#'
#' @examples
dml.lm.wrapper <- function(formula
                           , data
                           , first_stage_family = c('ols', 'rf', 'lasso', 'user-defined')
                           , predict_fun = NULL
                           , second_stage_family = c('mr', 'sr1', 'sr2')
                           , ...){

  # data.new = data
  f0 = Formula::Formula(formula)

  # y_vars (outcomes)
  y_vars = attr(terms(f0, lhs=T, rhs=0), 'term.labels')
  if(length(y_vars != 1)){
    stop('dml.lm.wrapper: one y_var at a time')
  }

  # x_vars (covariates)
  x_vars = attr(terms(f0, lhs=F, rhs=3), 'term.labels')
  if(any(!(x_vars %in% colnames(data)))){
    stop('dml.lm.wrapper: x_vars not in data')
  }

  # d_vars (treatments)
  d_vars = attr(terms(f0, lhs=F, rhs=2), 'term.labels')
  if(any(!(d_vars %in% colnames(data)))){
    stop('dml.lm.wrapper: d_vars not in data')
  }


  rhs.term.factors = attr(terms(f0, lhs=F, rhs=1), 'factors')
  treat.term.factors = rhs.term.factors[(rownames(rhs.term.factors) %in% d_vars), ]

  # h_vars (treatments and interactions)
  h_vars = data.frame(term=colnames(treat.term.factors))
  if(any(apply(treat.term.factors != 0, MARGIN=2, FUN=sum) != 1)){
    stop('dml.lm.wrapper: step 2 term missing a d_var or multiple d_vars')
  }
  h_vars$d = rownames(treat.term.factors)[apply(treat.term.factors != 0, MARGIN=2, FUN=which)]
  h_vars$fx = str_replace_all(h_vars$term, h_vars$d, '(Intercept)')
  h_vars$fxd.name = gsub('[:*]', '_', h_vars$term)

  # add transformations of covariates to data
  x0 = gsub('^:', '', gsub(':?\\(Intercept\\)', '', gsub('\\(Intercept\\):', '', h_vars$fx)))
  x1 = unique(x0)[unique(x0) != '']
  y0 = attr(terms(f0, lhs=T, rhs=0), 'term.labels')
  f1 = formula(paste(paste(y0, collapse=' + '), '~', paste(x1, collapse=' + ')))
  print(paste('dml.lm.wrapper: attempting to derive [', paste(x1, collapse=', '),
              '] from [', paste(x_vars, collapse=', '), ']', sep=''))
  data.new = cbind(data, model.matrix(f1, data=data[, x_vars]))

  ret <- dml.lm(data.new, y_var=y_vars, d_vars = NULL, x_vars = x_vars, h_vars = h_vars,
                first_stage_family=first_stage_family,
                predict_fun=predict_fun,
                second_stage_family=second_stage_family,
                ... = ...)
  return(ret)
}
