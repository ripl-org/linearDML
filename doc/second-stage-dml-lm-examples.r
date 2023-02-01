
data(iris)

#
# Using d_vars, no covariate interaction terms
#

# dml with OLS, RF prediction
for(fam in c('ols', 'rf')){
  model.dml.0 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width', 'Petal.Width'),
                       d_vars='Petal.Length', first_stage_family=fam)
  print(summary(model.dml.0$model)$coefficients)
}

# linear model, for comparison
model.lm.0 = lm(Sepal.Length ~ Petal.Length + Sepal.Width + Petal.Width, data=iris)
print(summary(model.lm.0)$coefficients)


#
# Using h_vars, yes covariate interaction terms
#
iris$sep.wid.2 = with(iris, Sepal.Width ^ 2)
iris$pet.len.sep.wid.2 = with(iris, Petal.Length * sep.wid.2)
iris$const = 1

h_vars = data.frame(d=c('Petal.Length', 'Petal.Length', 'Petal.Width'),
                    fx=c('const', 'sep.wid.2', 'const'),
                    fxd.name=c('Petal.Length', 'pet.len.sep.wid.2', 'Petal.Width'))

# dml with OLS, RF prediction
for(fam in c('ols', 'rf')){
  model.dml.1 = dml.lm(iris, y_var='Sepal.Length', x_vars=c('Sepal.Width'),
                       h_vars=h_vars, first_stage_family=fam)
  print(summary(model.dml.1$model)$coefficients)
}

# linear model, for comparison
model.lm.0 = lm(Sepal.Length ~ Petal.Length + Sepal.Width +
                  Petal.Width + pet.len.sep.wid.2, data=iris)
print(summary(model.lm.0)$coefficients)
