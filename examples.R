library(causaldata)
library(dplyr)
library(tidyr)
library(linearDML)

df <- nsw_mixtape %>%
  as.data.frame() %>%
  dplyr::mutate(treat =  as.logical(treat)
         , not_treated = !treat
         , nodegree_treated = treat*nodegree
         , degree_treated = treat*!nodegree)

#use OLS with multiple residualization
dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat')
       , first_stage_family = 'ols'
       , second_stage_family = 'mr')

#use OLS with single residualization
dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat', 'degree_treated')
       , first_stage_family = 'ols'
       , second_stage_family = 'sr1')

#use OLS with single residualization2
dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat', 'degree_treated')
       , first_stage_family = 'ols'
       , second_stage_family = 'sr2')

#use random forest with multiple residualization
dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat')
       , first_stage_family = 'rf'
       , second_stage_family = 'mr')

#use random forest with two hyperparameters set
dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat')
       , first_stage_family = 'rf'
       , second_stage_family = 'mr'
       , mtry = 3
       , trees = 100)

#use random forest with all hyperparameters set
dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat')
       , first_stage_family = 'rf'
       , second_stage_family = 'mr'
       , mtry = 3
       , trees = 100
       , max_depth = 2
       , min_n = 2)


#use a user defined prediction function
nnet_predict_fun <- function(data, y_var, x_vars){
  y <- data %>% pull(y_var)

  formula <- as.formula(paste(y_var, paste(x_vars, collapse=" + "), sep=" ~ "))

  model <- neuralnet::neuralnet(formula, data)
  y_hat <- predict(model, newdata = data)

  resids <- y - y_hat

  list(resids = resids)
}


dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('treat')
       , first_stage_family = 'user-defined'
       , predict_fun = nnet_predict_fun)
