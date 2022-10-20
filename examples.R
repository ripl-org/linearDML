library(causaldata)
library(dplyr)
library(tidyr)

devtools::load_all()

df <- nsw_mixtape %>%
  as.data.frame() %>%
  mutate(treat =  as.logical(treat)
         , not_treated = !treat
         , nodegree_treated = treat*nodegree
         , degree_treated = treat*!nodegree)

dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('not_treated')
       , first_stage_family = 'ols'
       , second_stage_family = 'mr')

dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('not_treated', 'degree_treated')
       , first_stage_family = 'ols'
       , second_stage_family = 'sr1')

dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('not_treated', 'degree_treated')
       , first_stage_family = 'ols'
       , second_stage_family = 'sr2')


dml.lm(data = df
       , x_vars = c('age', 'black', 'hisp', 'marr', 're74', 're75')
       , y_var = 're78'
       , d_vars = c('not_treated')
       , first_stage_family = 'rf'
       , second_stage_family = 'mr')
