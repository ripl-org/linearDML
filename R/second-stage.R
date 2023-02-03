
#' dml.lm(): Fitting double/debiased machine learning (DML) models
#'
#' Estimates treatment effects by regressing outcome residuals on treatment
#' residuals, where residuals come from predicting outcomes and treatments with
#' covariates.
#'
#' @importFrom Rdpack reprompt
#'
#' @param data A data frame.
#' @param y_var The name of the outcome variable in `data`. A character string.
#' @param d_vars The names of the treatment variables in `data`. A vector of character strings.
#' Either `d_vars` or `h_vars` must be specified.
#' If `h_vars` is specified, `d_vars` will be ignored.
#' @param x_vars The names of the covariates in `data`. A vector of character strings.
#' @param h_vars A data frame containing the following fields: `d`, `fx`, and `fxd.name`.
#' Each row represents a treatment variable interacted with a transformation of covariates.
#' `d` is the name of the corresponding treatment variable in `data`.
#' `fx` is the name of the covariate transformation in `data`.
#' If this covariate transformation is not actually a transformation of a
#' variable or set of variables in `x_vars`, then this function will not produce
#' consistent estimates.
#' `fxd.name` will be the name of the resulting interaction term (must be unique).
#' Either `d_vars` or `h_vars` must be specified.
#' If `h_vars` is specified, `d_vars` will be ignored.
#' @param first_stage_family The type of prediction done in the first stage.
#' Must be `'ols'` (Ordinary Least Squares), `'rf'` (Random Forest), or `'user-defined'`. If `'user-defined'`, then
#' `predict_fun` must also be specified. Note that OLS prediction will not
#' produce identical estimates to a long OLS regression including treatments
#' and covariates because DML uses out-of-sample prediction.
#' @param predict_fun The prediction function used in the first stage.
#' It must take as parameters `data` (a data frame), `y_var` (name of outcome or treatment),
#' and `x_vars` (names of covariates). It must return a list with a field
#' `resids` that contains the residuals from predicting `y_var` using `x_vars`
#' as a numeric vector of length `nrow(data)`.
#' If `first_stage_family` == `'user-defined'`, then `predict_fun` must be specified.
#' Otherwise, `predict_fun` will be ignored.
#' @param second_stage_family Indicates how treatments are residualized.
#' Can be `'mr'` (Multiple Residualization, the standard DML method),
#' `'sr1'` (Single Residualization 1), or `'sr2'` (Single Residualization 2).
# @param use_default_hyper placeholder
#'
#' @return `dml.lm` returns an list containing the following components: \cr
#' \tabular{ll}{
#'  `model`       \tab The `lm` object used to fit the second stage. The coefficients from this model are second-stage DML estimates. \cr
#'  `y_model`     \tab The object used to residualize the outcome variable `y_var` in the first stage. \cr
#'  `d_model`     \tab The object used to residualize the treatment variables `d_vars` in the first stage.
#' }
#'
#' @details `dml.lm` estimates Double Machine Learning (DML) models
#' according to \insertCite{chernozhukov_doubledebiased_2018;textual}{riplDML}.
#' In that framework, both treatment and outcome variables are residualized by
#' a prediction function (the first stage), and then the outcome residuals are regressed on the
#' treatment residuals (the second stage). \cr
#' This procedure is able to control for large numbers of covariates (and
#' covariate transformations) by incorporating variable selection (e.g.
#' through LASSO or Random Forest) in the first stage. The estimates are
#' consistent because both the treatments and the covariates have been
#' residualized. \cr
#' We allow for treatment effect heterogeneity according to the procedure
#' summarized in \insertCite{semenova2017estimation;textual}{riplDML}. When defined, the
#' parameter `h_vars` contains the set of treatments including covariate
#' transformation interaction terms. \cr
#' This function returns an `lm` object with default homoskedastic standard
#' errors. That object can be passed into a `vcov`-style function to get
#' heteroskedasticity-robust and/or cluster-robust standard errors. The
#' errors from the second stage correspond to the score function estimators
#' from \insertCite{chernozhukov_doubledebiased_2018;textual}{riplDML}, and the
#' residualized outcome and treatments correspond to the components of the
#' linear score.
#'
#' @references{
#'  \insertRef{chernozhukov_doubledebiased_2018}{riplDML} \cr
#'  \insertRef{semenova2017estimation}{riplDML}
#' }
#'
#' @export
#'
#' @examples
#'
#' data(iris)
#'
#' #
#' # Using d_vars, no covariate interaction terms
#' #
#'
#' # dml with OLS, RF prediction
#' for(fam in c('ols', 'rf')){
#'   model.dml.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
#'                        d_vars='Petal.Length', first_stage_family=fam)
#'   print(summary(model.dml.0$model)$coefficients)
#' }
#'
#' # linear model, for comparison
#' model.lm.0 = lm(Sepal.Length ~ Petal.Length + Sepal.Width + Petal.Width, data=iris)
#' print(summary(model.lm.0)$coefficients)
#'
#'
#' #
#' # Using h_vars, yes covariate interaction terms
#' #
#' iris$sep.wid.2 = with(iris, Sepal.Width ^ 2)
#' iris$pet.len.sep.wid.2 = with(iris, Petal.Length * sep.wid.2)
#' iris$const = 1
#'
#' h_vars = data.frame(d=c('Petal.Length', 'Petal.Length', 'Petal.Width'),
#'                     fx=c('const', 'sep.wid.2', 'const'),
#'                     fxd.name=c('Petal.Length', 'pet.len.sep.wid.2', 'Petal.Width'))
#'
#' # dml with OLS, RF prediction
#' for(fam in c('ols', 'rf')){
#'   model.dml.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width'),
#'                        h_vars=h_vars, first_stage_family=fam)
#'   print(summary(model.dml.1$model)$coefficients)
#' }
#'
#' # linear model, for comparison
#' model.lm.0 = lm(Sepal.Length ~ Petal.Length + Sepal.Width +
#'                   Petal.Width + pet.len.sep.wid.2, data=iris)
#' print(summary(model.lm.0)$coefficients)
dml.lm <- function(data
                   , y_var
                   , x_vars
                   , d_vars = NULL
                   , h_vars = NULL
                   , first_stage_family = c('ols', 'rf', 'user-defined')
                   , predict_fun = NULL
                   , second_stage_family = c('mr', 'sr1', 'sr2')
                   , ...){
  out = list()

  ## evaluate choices
  first_stage_family <- match.arg(first_stage_family)
  second_stage_family <- match.arg(second_stage_family)

  if(is.null(d_vars)){
    if(second_stage_family %in% c('sr1', 'sr2')){
      stop('dml.lm: if second_stage_family is in (sr1, sr2), d_vars cannot be NULL')
    }
    if(is.null(h_vars)){
      stop('dml.lm: if d_vars and h_vars cannot both be NULL')
    }
  }

  if(length(y_var) != 1){
    stop('dml.lm: length(y_var) != 1')
  }

  if(!is.null(d_vars)){
    factor_d_names <- names(dplyr::select(data, d_vars))[ sapply(dplyr::select(data, d_vars) , is.factor)]
    non_factor_d_names <- names(dplyr::select(data, d_vars))[ sapply(dplyr::select(data, d_vars), function(x) !is.factor(x))]

    if(length(factor_d_names) > 0){
      factor_names_str <- paste(factor_d_names, sep = ',')
      message(paste0('Converting factors ', factor_names_str, ' to dummy columns.'))

      old_colnames <- colnames(data)
      data <- data %>% fastDummies::dummy_cols(factor_d_names)
      new_factor_cols <- setdiff(colnames(data), old_colnames)

      #reform d_vars to use the dummy column names
      d_vars <- c(non_factor_d_names, new_factor_cols)

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

      tmp1 = duplicated(h_vars$fxd.name)
      if(any(tmp1)){
        stop(paste('h_vars$fxd.name must be unique;',
                   h_vars$fxd.name[tmp1][1], 'is duplicated'))
      }

      # Predict treatments
      d_models <- lapply(d_vars, function(d_var) predict.dml(data
                                                             , y_var = d_var
                                                             , x_vars = x_vars
                                                             , family = first_stage_family
                                                             , ...))

      d_resids <- as.data.frame(do.call(cbind, d_models %>% lapply(purrr::pluck, 'resids')))
      names(d_resids) <- d_vars
      out$d_model <- d_models

            #
      # h_vars should be a data frame with three columns:
      # d - the treatment variable
      # fx - the (transformation of the) covariates
      # fxd.name - the name of this combined variable in the second stage
      reg.data <- data.frame(y_resids)
      # reg.data[, h_vars$fxd.name] <- data[, h_vars$fx] * purrr::pluck(d_resids, h_vars$d)
      reg.data[, h_vars$fxd.name] <- data[, h_vars$fx] * d_resids[, h_vars$d]
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
#' @param formula placeholder
#' @param data placeholder
#' @param first_stage_family placeholder
#' @param second_stage_family placeholder
#' @param use_default_hyper placeholder
#'
#' @return placeholder
# @export
#'
#' @examples placeholder
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
