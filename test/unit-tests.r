


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

  print('test_dml.lm_first_stage_family test completed')
}


#======================================================
# RUN ALL TESTS
#======================================================

test_dml.lm_y_var()
test_dml.lm_x_vars()
test_dml.lm_d_vars()
test_dml.lm_first_stage_family()
test_dml.lm_second_stage_family()








