START_TIME = Sys.time()

library(dplyr)
library(reshape)
library(stringr)
library(DoubleML)
library(mlr3)
library(mlr3learners)

# Description:
# This file demonstrates how dml.lm() with h_vars specified leverages a single
# treatment residualization across multiple treatment heterogeneity terms
# to reduce the computational burden of estimation.

# Load the package
devtools::load_all()
# Find the current script
script.name = basename(sys.frame(1)$ofile)
script.dir = dirname(sys.frame(1)$ofile)

lgr::get_logger("mlr3")$set_threshold("warn")

# Simulation parameters
NTREES = 50
MIN_N = 10
NSIM = 10
NOBS = 5000

# Generate a matrix with simulation and observation IDs
sim.0 = expand.grid(id=1:NOBS, sim=1:NSIM)

# Set a random seed
set.seed(20230216)

sim.0$one = 1
sim.0$c0 = 1

# functions for generating a few non-linear transformations
distr.0 = list()
distr.0[[1]] = function(x){ return(c(scale(qnorm(x)))) }
distr.0[[2]] = function(x){ return(c(scale(qexp(x)))) }
distr.0[[3]] = function(x){ return(c(scale(qbeta(x, 0.1, 0.2)))) }
distr.0[[4]] = function(x){ return(c(scale(qbinom(x, size=20, prob=0.08)))) }

# structure many separate simulations
sim.0$big.id = 1:nrow(sim.0)
sim.1 = sim.0 %>%
  group_by(sim) %>%
  summarize(obs=sum(one),
            ind.lo=min(big.id),
            ind.hi=max(big.id),
            .groups='keep')
sim.1 = as.data.frame(sim.1)
rownames(sim.1) = paste('tag', sim.1$sim, 'end', sep='.')

# primitives
for(i in 1:20){
  sim.1[, paste('delta1', i, sep='.')] = rnorm(nrow(sim.1))
  del1 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), paste('delta1', i, sep='.')]
  sim.0[, paste('delta1', i, sep='.')] = del1

  sim.1[, paste('delta2', i, sep='.')] = rnorm(nrow(sim.1))
  del2 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), paste('delta2', i, sep='.')]
  sim.0[, paste('delta2', i, sep='.')] = del2

  sim.1[, paste('gamma', i, sep='.')] = rnorm(nrow(sim.1))
  gam = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), paste('gamma', i, sep='.')]
  sim.0[, paste('gamma', i, sep='.')] = gam

  tmp1 = runif(nrow(sim.0))
  sim.0[, paste('v1', i, sep='.')] = distr.0[[(i %% 4) + 1]](tmp1) * del1
  sim.0[, paste('v2', i, sep='.')] = distr.0[[(i %% 4) + 1]](tmp1) * del2
  sim.0[, paste('w', i, sep='')] = ((-1)^(i %% 3)) * distr.0[[((i + 1) %% 4) + 1]](tmp1)

  sim.0[, paste('u', i, sep='')] = distr.0[[(i %% 4) + 1]](tmp1) * gam
}

sim.0$x1.0 = apply(sim.0[, paste('v1', 1:10, sep='.')], MARGIN=1, FUN=mean)
sim.0$x2.0 = apply(sim.0[, paste('v2', 6:15, sep='.')], MARGIN=1, FUN=mean)

# binary treatment
sig = sd(sim.0$x1.0)
sim.0$eta1 = rnorm(nrow(sim.0), mean=sig, sd=sig * 0.9)
sim.0$x1 = as.numeric((sim.0$x1.0 + sim.0$eta1) < 0)

# continuous treatment
sim.0$eta2 = rnorm(nrow(sim.0), mean=sig, sd=sig * 0.9)
sim.0$x2 = sim.0$x2.0 + sim.0$eta2

# diagnostic
look.1 = sim.0[sim.0$sim == 5, ]
for(i in 1:20){
  plot(look.1[, paste('w', i, sep='')], look.1[, paste('u', i, sep='')], main=i)
}


#
# construct outcomes
#

# sum up covariate effects
sim.0$u.all = apply(sim.0[, paste('u', 1:20, sep='')], MARGIN=1, FUN=sum)
# exogenous errors
sim.0$e1 = rnorm(nrow(sim.0)) + (rbinom(nrow(sim.0), size=2, prob=0.5) * 4)
# constant
gamma.0 = 3

# treatment effects - main
sim.1$alpha1 = rnorm(nrow(sim.1))
sim.0$alpha1 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'alpha1']

sim.1$alpha2 = rnorm(nrow(sim.1))
sim.0$alpha2 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'alpha2']

# treatment effects - heterogeneous
sim.0$het.all = 0
for(i in 1:20){
  for(j in 1:2){
    pn = paste('beta', j, '.', i, sep='')
    sim.1[, pn] = rnorm(nrow(sim.1))
    sim.0[, pn] = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), pn]

    tn = paste('x', j, sep='')
    vn = paste('w', i, sep='')
    sim.0$het.all = sim.0$het.all + (sim.0[, pn] * sim.0[, tn] * sim.0[, vn])

    hn = paste(tn, vn, sep='_')
    sim.0[, hn] = sim.0[, tn] * sim.0[, vn]
  }
}

# 1. no heterogeneity, only x1
sim.0$y1 = with(sim.0, gamma.0 + u.all + (alpha1 * x1) + e1)

# 2. yes heterogeneity, only x1
sim.0$y2 = with(sim.0, gamma.0 + u.all + (alpha1 * x1) +
                  (beta1.10 * x1 * w10) + (beta1.20 * x1 * w20) + e1)

# 3. no heterogeneity, only x2
sim.0$y3 = with(sim.0, gamma.0 + u.all + (alpha2 * x2) + e1)

# 4. yes heterogeneity, only x2
sim.0$y4 = with(sim.0, gamma.0 + u.all + (alpha2 * x2) +
                  (beta2.10 * x2 * w10) + (beta2.20 * x2 * w20) + e1)

# 5. yes heterogeneity, all treatments and covariates
sim.0$y5 = with(sim.0, gamma.0 + u.all + (alpha1 * x1) + (alpha2 * x2) +
                  het.all + e1)

keeps = c('sim', 'c0', paste('y', 1:5, sep=''), 'x1', 'x2', paste('w', 1:20, sep=''),
          paste('x1_w', 1:20, sep=''), paste('x2_w', 1:20, sep=''))
sim.2 = sim.0[, keeps]

#===========================================================
# Estimation Start
#===========================================================

# Specifying h_vars for each of the 5 models
h.vars.1 = data.frame(d=c('x1'),
                      fx=c('c0'),
                      fxd.name=c('x1'))
h.vars.2 = data.frame(d=c('x1', 'x1', 'x1'),
                      fx=c('c0', 'w10', 'w20'),
                      fxd.name=c('x1', 'x1_w10', 'x1_w20'))
h.vars.3 = data.frame(d=c('x2'),
                      fx=c('c0'),
                      fxd.name=c('x2'))
h.vars.4 = data.frame(d=c('x2', 'x2', 'x2'),
                      fx=c('c0', 'w10', 'w20'),
                      fxd.name=c('x2', 'x2_w10', 'x2_w20'))
h.vars.5 = data.frame(d=c(rep('x1', 21), rep('x2', 21)),
                      fx=c('c0', paste('w', 1:20, sep=''), 'c0', paste('w', 1:20, sep='')),
                      fxd.name=c('x1', paste('x1_w', 1:20, sep=''), 'x2', paste('x2_w', 1:20, sep='')))
w.all = paste('w', 1:20, sep='')


#
# Estimation of one simulation at a time
# using linearDML
#
coeff.10 = NULL
for(i in 1:NSIM){
  print(paste(Sys.time(), i, sep=' : '))

  sim.3 = sim.2[sim.2$sim == sim.1$sim[i], ]

  for(j in 1:5){
    h.vars = get(paste('h.vars', j, sep='.'))
    y.curr = paste('y', j, sep='')
    d.vars = h.vars$fxd.name

    #
    # DML with OLS prediction
    #
    TIME_00 = Sys.time()
    model.0 = (dml.lm(sim.3, y.curr, w.all, h_vars=h.vars,
                                      first_stage_family = 'ols', second_stage_family = 'mr'))
    coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
    colnames(coeff.0) = c('b', 'se', 't', 'p')
    coeff.0$x = rownames(coeff.0)
    coeff.0$method = 'ols'
    coeff.0$y = j
    coeff.0$sim = i
    coeff.10 = rbind(coeff.10, coeff.0)
    TIME_01 = Sys.time()
    sim.1[i, paste('time.ols', j, sep='.')] = TIME_01 - TIME_00

    #
    # DML with RF prediction
    #
    TIME_00 = Sys.time()
    model.0 = suppressMessages(dml.lm(sim.3, y.curr, w.all, h_vars=h.vars,
                                      first_stage_family = 'rf', second_stage_family = 'mr',
                                      trees = NTREES, min_n=MIN_N, mtry=10, max_depth=5,
                                      use_default_hyper = F))
    coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
    colnames(coeff.0) = c('b', 'se', 't', 'p')
    coeff.0$x = rownames(coeff.0)
    coeff.0$method = 'rf'
    coeff.0$y = j
    coeff.0$sim = i
    coeff.10 = rbind(coeff.10, coeff.0)
    TIME_01 = Sys.time()
    sim.1[i, paste('time.rf', j, sep='.')] = TIME_01 - TIME_00

    if(T){
      #
      # DML with RF prediction, not calculated heterogeneously
      #
      TIME_00 = Sys.time()
      model.0 = suppressMessages(dml.lm(sim.3, y.curr, w.all, d_vars=d.vars,
                                        first_stage_family = 'rf', second_stage_family = 'mr',
                                        trees = NTREES, min_n=MIN_N, mtry=10, max_depth=5,
                                        use_default_hyper = F))
      coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
      colnames(coeff.0) = c('b', 'se', 't', 'p')
      coeff.0$x = rownames(coeff.0)
      coeff.0$method = 'rf2'
      coeff.0$y = j
      coeff.0$sim = i
      coeff.10 = rbind(coeff.10, coeff.0)
      TIME_01 = Sys.time()
      sim.1[i, paste('time.rf2', j, sep='.')] = TIME_01 - TIME_00


      # DML from DoubleML

      TIME_00 = Sys.time()
      learner = lrn("regr.ranger", num.trees = NTREES, min.node.size=MIN_N, mtry=10, max.depth=5)
      ml_l = learner$clone()
      ml_m = learner$clone()
      obj_dml_data = suppressMessages(DoubleMLData$new(sim.3[, c(y.curr, d.vars, w.all)],
                                                       y_col=y.curr, d_cols=d.vars))
      dml_plr_obj = suppressMessages(DoubleMLPLR$new(obj_dml_data, ml_l, ml_m, n_folds=2))
      suppressMessages(dml_plr_obj$fit())

      invisible(capture.output(coeff.0 <- as.data.frame(as.matrix(dml_plr_obj$summary()))))
      colnames(coeff.0) = c('b', 'se', 't', 'p')
      coeff.0$x = rownames(coeff.0)
      coeff.0$method = 'dub'
      coeff.0$y = j
      coeff.0$sim = i
      coeff.10 = rbind(coeff.10, coeff.0)
      TIME_01 = Sys.time()
      sim.1[i, paste('time.dub', j, sep='.')] = TIME_01 - TIME_00
    }


  }
}

# Examine main effects and heterogeneity terms
params.0 = names(sim.1)[grepl('((alpha)|(beta))', names(sim.1))]
params.1 = as.data.frame(melt(sim.1, id.vars=c('sim'),
                              measure.vars=params.0))

tmp1 = grepl('^alpha[12]$', params.1$variable)
params.1[tmp1, 'x'] = paste('x', substr(params.1[tmp1, 'variable'], 6, 6), sep='')
tmp1 = grepl('^beta[12].[0-9]{1,2}$', params.1$variable)
params.1[tmp1, 'x'] = paste('x', substr(params.1[tmp1, 'variable'], 5, 5),
                            '_w', substr(params.1[tmp1, 'variable'], 7, 8), sep='')
params.1$param = params.1$value

keeps = c('sim', 'x', 'param')
params.2 = params.1[, keeps]
coeff.11 = merge(coeff.10, params.2, by=c('sim', 'x'), all.x=T, all.y=F)
coeff.12 = coeff.11[with(coeff.11, order(sim, y, method)), ]

coeff.12$one = 1

# get mean absolute distance and root mean squared error across all simulations
# for each model.
coeff.13 = coeff.12 %>%
  filter(sim <= 5) %>%
  group_by(y, x, mm=method) %>%
  summarize(obs=sum(one),
            mad=mean(abs(b - param)),
            rmse=sqrt(mean((b - param)^2)),
            .groups='keep')
coeff.13 = as.data.frame(coeff.13)
coeff.14 = as.data.frame(cast(coeff.13, y + x + obs ~ mm, value='rmse'))
coeff.14$w = as.numeric(ifelse(grepl('w', coeff.14$x), gsub('x[0-9]_w', '', coeff.14$x), '0'))
coeff.14$x0 = gsub('_w[0-9]{1,2}$', '', coeff.14$x)
coeff.15 = coeff.14[with(coeff.14, order(y, x0, w, obs)), ]

# Examine time elapsed in each method
time.0 = names(sim.1)[grepl('(time)', names(sim.1))]
time.1 = as.data.frame(melt(sim.1, id.vars=c('sim'),
                              measure.vars=time.0))
time.1$method = gsub('\\.[0-9]$', '', gsub('^time\\.', '', time.1$variable))
time.1$y = gsub('^[a-z0-9]{1,5}\\.', '', gsub('^time\\.', '', time.1$variable))
time.1$one = 1

time.2 = time.1 %>%
  group_by(y, mm=method) %>%
  summarize(obs=sum(one),
            time.sum=sum(value),
            .groups='keep')
time.2 = as.data.frame(time.2)

time.3 = as.data.frame(cast(time.2, y + obs ~ mm, value='time.sum'))

print(script.dir)
setwd(script.dir)

write.csv(coeff.10, 'test-hetcoefall.csv', row.names=F)
write.csv(coeff.15, 'test-hetcoef.csv', row.names=F)
write.csv(time.3, 'test-hettime.csv', row.names=F)


#===========================================================
# Estimation End
#===========================================================

END_TIME = Sys.time()
print(paste('Time elapsed in ', script.name, ':', sep=''))
print(END_TIME - START_TIME)



