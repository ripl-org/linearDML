library(linearDML)
library(causaldata)
library(dplyr)
library(tidyr)
library(neuralnet)

test_dml.lm_y_var <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')
  model.1 = dml.lm(iris, y_var='Sepal.Width', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')

  print('test_dml.lm_y_var test completed')
}



test_dml.lm_x_vars <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')
  model.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width', 'Petal.Length'),
                   d_vars='Petal.Length', first_stage_family='ols')

  print('test_dml.lm_x_vars test completed')
}


test_dml.lm_d_vars <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')
  model.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Species', first_stage_family='ols')



  print('test_dml.lm_d_vars test completed')
}

test_dml.lm_first_stage_family <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')
  model.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length')
  model.2 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='rf')

  print('test_dml.lm_first_stage_family test completed')
}

test_dml.lm_second_stage_family <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Species', second_stage_family='mr')
  model.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Species', second_stage_family='sr1')
  model.2 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Species', second_stage_family='sr2')
  model.3 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Species')
  model.4 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Species', first_stage_family='rf', second_stage_family='sr2')

  print('test_dml.lm_second_stage_family test completed')
}

test_dml.lm_h_vars <- function(){
  data(iris)

  iris$sep.wid.2 = with(iris, Sepal.Width ^ 2)
  iris$pet.len.sep.wid.2 = with(iris, Petal.Length * sep.wid.2)
  iris$const = 1

  h_vars = data.frame(d=c('Petal.Length', 'Petal.Length', 'Petal.Width'),
                      fx=c('const', 'sep.wid.2', 'const'),
                      fxd.name=c('Petal.Length', 'pet.len.sep.wid.2', 'Petal.Width'))
  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width'),
                   h_vars=h_vars, first_stage_family='ols')

  h_vars = data.frame(d=c('Petal.Length'),
                      fx=c('const'),
                      fxd.name=c('Petal.Length'))
  model.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   h_vars=h_vars, first_stage_family='ols')

  print('test_dml.lm_h_vars test completed')
}

test_dml.hyper_parameters <- function(){
  df <- nsw_mixtape %>%
    as.data.frame() %>%
    dplyr::mutate(treat =  as.logical(treat))


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
}

test_dml.user_function <- function(){

  df <- nsw_mixtape %>%
    as.data.frame() %>%
    dplyr::mutate(treat =  as.logical(treat))

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

}

#======================================================
# RUN ALL TESTS
#======================================================

test_dml.lm_y_var()
test_dml.lm_x_vars()
test_dml.lm_d_vars()
test_dml.lm_first_stage_family()
test_dml.lm_second_stage_family()
test_dml.lm_h_vars()
test_dml.hyper_parameters()
test_dml.user_function()




