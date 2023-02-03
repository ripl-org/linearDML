


test_dml.lm_y_var <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')
  model.1 = dml.lm(iris, y_var='Sepal.Width', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')

  print('test completed')
}



test_dml.lm_y_var <- function(){
  data(iris)

  model.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')
  model.1 = dml.lm(iris, y_var='Sepal.Width', x_vars=c('Sepal.Width', 'Petal.Width'),
                   d_vars='Petal.Length', first_stage_family='ols')

  print('test completed')
}






#======================================================
# RUN ALL TESTS
#======================================================

test_dml.lm_y_var()








