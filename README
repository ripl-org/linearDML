# linearDML

## Overview

`linearDML` implements double/debiased machine learning (DML) Chernozhukov et al. (2018) and adds a method for flexibly specifying interactions with treatment variables.

## Usage

`dml.lm` is the workhorse function of `linearDML`. Along with first and second stage DML estimates, `dml.lm` always returns an `lm` object with the average treatment effect

```r
# dml with RF predictions

model.dml = dml.lm(iris
  , y_var = 'Sepal.Length' #outcome
  , x_vars = c('Sepal.Width', 'Petal.Width') #co-variates
  , d_vars = 'Petal.Length' #treatment
  , first_stage_family = 'rf')

summary(model.dml$model)
```
When using `dml.lm` to specify variable interactions, the `h_vars` argument expects a dataframe with three columns

* `d` - the name of the treatment variable to interact
* `fx` - the name of the covariate to interact
* `fxd.name` - the name of a calculated columns that defines the variable interaction

```
# Using h_vars, with covariate interaction terms
#
iris$sep.wid.2 = with(iris, Sepal.Width ^ 2)
iris$pet.len.sep.wid.2 = with(iris, Petal.Length * sep.wid.2)
iris$const = 1


h_vars = data.frame(d = c('Petal.Length', 'Petal.Length', 'Petal.Width'),
                    fx = c('const', 'sep.wid.2', 'const'),
                    fxd.name = c('Petal.Length', 'pet.len.sep.wid.2', 'Petal.Width'))

model.dml.h = dml.lm(iris
  , y_var = 'Sepal.Length'
  , x_vars = c('Sepal.Width')
  , h_vars = h_vars
  , first_stage_family = 'ols')

summary(model.dml.h$model)

```

## Installation

You can install the latest version of `linearDML` from Github using the following command

```r
devtools::install_github('ripl-org/linearDML')
```


## License
[LICENSE](LICENSE)

Copyright 2023 Innovative Policy Lab d/b/a Research Improving People's Lives
("RIPL"), Providence, RI. All Rights Reserved.

## References
