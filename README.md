# linearDML

## Overview

`linearDML` implements double/debiased machine learning (DML) [Chernozhukov et al. (2018)](https://doi.org/10.1111/ectj.12097) for partially linear models and adds a method for flexibly specifying interactions with treatment variables.


## Usage

`dml.lm` is the workhorse function of `linearDML`. Along with first stage DML predictions, `dml.lm` returns an `lm` object containing estimates from the second stage.

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

```r
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

## References

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018),
Double/debiased machine learning for treatment and structural parameters. The Econometrics Journal, 21: C1-C68. doi:[10.1111/ectj.12097](https://doi.org/10.1111/ectj.12097).


## License
[LICENSE](LICENSE)

Copyright 2023 Innovative Policy Lab d/b/a Research Improving People's Lives
("RIPL"), Providence, RI. All Rights Reserved.
