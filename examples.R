library(causaldata)
library(dplyr)
library(tidyr)

devtools::load_all()

dml_first_stage_rf(data = nsw_mixtape
                   , x_vars = c(age, educ, black, hisp, marr, re74, re75)
                   , y_var = re78
                   , d_vars = c(treat, nodegree))


load_all()
dml_first_stage_lasso(data = nsw_mixtape
                   , x_vars = c(age, educ, black, hisp, marr, re74, re75)
                   , y_var = re78
                   , y_family = 'gaussian'
                   , d_vars = c(treat, nodegree)
                   , d_family = 'binomial')
