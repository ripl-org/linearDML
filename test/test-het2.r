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
NSIM = 100
NOBS = 5000
MTRY = 4

# Generate a matrix with simulation and observation IDs
sim.0 = expand.grid(id=1:NOBS, sim=1:NSIM)

# Set a random seed
set.seed(20230216)

sim.0$one = 1
sim.0$c0 = 1

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

# generate parameters for each simulation. these parameters are constant
# within each simulation but vary across simulations.
sim.1$delta = runif(nrow(sim.1)) + 1
sim.0$delta = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'delta']
sim.1$gamma = runif(nrow(sim.1)) + 1
sim.0$gamma = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'gamma']

sim.0$wallmid = 0
for(i in 1:5){
  tmp1 = runif(nrow(sim.0))
  sim.0[, paste('w', i, sep='')] = tmp1
  sim.0[, paste('w', i, 'mid', sep='')] = as.numeric(((1/3) < tmp1) & (tmp1 < (2/3)))
  sim.0$wallmid = sim.0$wallmid + sim.0[, paste('w', i, 'mid', sep='')]
}

# Generate covariate transformations that are not correlated with any of the
# covariates themselves
sim.0$vcorners = with(sim.0, as.numeric((w1mid + w2mid) == 0))
sim.0$vedges = with(sim.0, as.numeric((w1mid + w2mid) == 1))
sim.0$vskel = with(sim.0, as.numeric((w1mid + w2mid) <= 1))
look.0 = cor(sim.0[, c('w1', 'w2', 'w3', 'w4', 'w5', 'vcorners', 'vedges', 'vskel')])

# Latent treatments
sim.0$x1.0 = sim.0$vcorners * sim.0$delta
sim.0$x2.0 = sim.0$vcorners * sim.0$delta

# binary treatment
sig = sd(sim.0$x1.0)
sim.0$eta1 = rnorm(nrow(sim.0), mean=0, sd=1)
sim.0$x1 = as.numeric((sim.0$x1.0 + sim.0$eta1) < 0)

# continuous treatment
sim.0$eta2 = rnorm(nrow(sim.0), mean=0, sd=1)
sim.0$x2 = sim.0$x2.0 + sim.0$eta2

#
# construct outcomes
#

sim.0$e1 = rnorm(nrow(sim.0)) + (rbinom(nrow(sim.0), size=2, prob=0.5) * 4)
gamma.0 = 3

# treatment effects - main
sim.1$alpha1 = rnorm(nrow(sim.1))
sim.0$alpha1 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'alpha1']

sim.1$alpha2 = rnorm(nrow(sim.1))
sim.0$alpha2 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'alpha2']

# treatment effects - heterogeneous
sim.0$het.all.1 = 0
sim.0$het.all.2 = 0
for(i in 1:5){
  for(j in 1:2){
    pn = paste('beta', j, '.', i, sep='')
    sim.1[, pn] = rnorm(nrow(sim.1))
    sim.0[, pn] = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), pn]

    tn = paste('x', j, sep='')
    vn = paste('w', i, sep='')
    hn = paste(tn, vn, sep='_')
    sim.0[, hn] = sim.0[, tn] * sim.0[, vn]

    sim.0[, paste('het.all', j, sep='.')] =
      sim.0[, paste('het.all', j, sep='.')] +
      sim.0[, pn] * sim.0[, hn]
  }
}

# 1. no heterogeneity, only x1
sim.0$y1 = with(sim.0, gamma.0 + (alpha1 * x1) + (gamma * vcorners) + e1)

# 2. yes heterogeneity, only x1
sim.0$y2 = with(sim.0, gamma.0 + (alpha1 * x1) + (gamma * vcorners) +
                  het.all.1 + e1)

# 3. no heterogeneity, only x2
sim.0$y3 = with(sim.0, gamma.0 + (alpha2 * x2) + (gamma * vcorners) + e1)

# 4. yes heterogeneity, only x2
sim.0$y4 = with(sim.0, gamma.0 + (alpha2 * x2) + (gamma * vcorners) +
                  het.all.2 + e1)

keeps = c('sim', 'c0', paste('y', 1:4, sep=''), 'x1', 'x2', paste('w', 1:5, sep=''),
          paste('x1_w', 1:5, sep=''), paste('x2_w', 1:5, sep=''))
sim.2 = sim.0[, keeps]

#===========================================================
# Estimation Start
#===========================================================

h.vars.1 = data.frame(d=c('x1'),
                      fx=c('c0'),
                      fxd.name=c('x1'))
h.vars.2 = data.frame(d=c(rep('x1', 6)),
                      fx=c('c0', paste('w', 1:5, sep='')),
                      fxd.name=c('x1', paste('x1_w', 1:5, sep='')))
h.vars.3 = data.frame(d=c('x2'),
                      fx=c('c0'),
                      fxd.name=c('x2'))
h.vars.4 = data.frame(d=c(rep('x2', 6)),
                      fx=c('c0', paste('w', 1:5, sep='')),
                      fxd.name=c('x2', paste('x2_w', 1:5, sep='')))
w.all = paste('w', 1:5, sep='')


#
# Estimation of one simulation at a time
# using linearDML
#
coeff.10 = NULL
for(i in 1:NSIM){
  if((i %% 10) == 1){
    print(paste(Sys.time(), i, sep=' : '))
  }

  sim.3 = sim.2[sim.2$sim == sim.1$sim[i], ]

  for(j in 1:4){
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
                                      trees = NTREES, min_n=MIN_N, mtry=MTRY, max_depth=5,
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
                                        trees = NTREES, min_n=MIN_N, mtry=MTRY, max_depth=5,
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



      #
      # DML from DoubleML
      #
      TIME_00 = Sys.time()
      learner = lrn("regr.ranger", num.trees = NTREES, min.node.size=MIN_N, mtry=MTRY, max.depth=5)
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
            bias=mean(b - param),
            rmse=sqrt(mean((b - param)^2)),
            .groups='keep')
coeff.13 = as.data.frame(coeff.13)
coeff.14 = as.data.frame(melt(coeff.13, id.vars=c('y', 'x', 'mm'),
                              measure.vars=c('mad', 'bias', 'rmse')))
coeff.15 = as.data.frame(cast(coeff.14, y + x ~ variable + mm, value='value'))

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

write.csv(coeff.10, 'test-het2coefall.csv', row.names=F)
write.csv(coeff.15, 'test-het2coef.csv', row.names=F)
write.csv(time.3, 'test-het2time.csv', row.names=F)


#===========================================================
# Estimation End
#===========================================================

END_TIME = Sys.time()
print(paste('Time elapsed in ', script.name, ':', sep=''))
print(END_TIME - START_TIME)



