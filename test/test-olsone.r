START_TIME = Sys.time()

library(dplyr)
library(reshape)
library(DoubleML)
library(mlr3)
library(mlr3learners)

# Description:
# This file replicates the supposed functionality of linearDML::dml.lm()
# when the prediction function is OLS or RF, running regressions one by one
# using lm() and DoubleML. Then it checks those estimates
# against the actual linearDML::dml.lm() function.

# Load the package
devtools::load_all()
# Find the current script
script.name = basename(sys.frame(1)$ofile)
script.dir = dirname(sys.frame(1)$ofile)

lgr::get_logger("mlr3")$set_threshold("warn")

# Simulation parameters
NTREES = 100
MIN_N = 2
NSIM = 50
NOBS = 500

# Generate a matrix with simulation and observation IDs
sim.0 = expand.grid(id=1:NOBS, sim=1:NSIM)

# Set a random seed
set.seed(20220822)

sim.0$one = 1
sim.0$c0 = 1

# primitives - all gaussian
sim.0$w1 = rnorm(nrow(sim.0)) + 1
sim.0$w2 = rnorm(nrow(sim.0)) + 1
sim.0$w3 = rnorm(nrow(sim.0))
sim.0$eta1 = rnorm(nrow(sim.0))
sim.0$eta2 = rnorm(nrow(sim.0))
sim.0$epsilon = rnorm(nrow(sim.0))

# data generating process for treatment variables, x's.
# importantly, the covariance matrix between the x's and
# the w's is full rank. the data generating process for
# x1 is not the same as that for x2 up to a constant.
g11 = 3
g12 = 3
g13 = 1
g21 = -3
g22 = 3
g23 = 1
sim.0$x1.star = with(sim.0, (g11 * w1) + (g12 * w2) + (g13 * eta1))
sim.0$x2.star = with(sim.0, (g21 * w1) + (g22 * w2) + (g23 * eta1))
sim.0$x1 = with(sim.0, (w3 >= 0) * (x1.star >= 0))
sim.0$x2 = with(sim.0, (w3 < 0) * (x2.star >= 0))

sim.0$xbar = with(sim.0, x1 + x2)

# As stated above, the covariance matrix between the x's and
# the w's is full rank. These three lines allow us to verify
# this by examining Delta_0.
W = as.matrix(sim.0[, c('c0', 'w1', 'w2')])
X = as.matrix(sim.0[, c('x1', 'x2')])
Delta_0 = solve(t(W) %*% W) %*% (t(W) %*% X)

# structure many separate simulations
sim.0$big.id = 1:nrow(sim.0)
sim.1 = sim.0 %>%
  group_by(sim) %>%
  summarize(obs=sum(one),
            x1=sum(x1),
            x2=sum(x2),
            xbar=sum(xbar),
            ind.lo=min(big.id),
            ind.hi=max(big.id),
            .groups='keep')
sim.1 = as.data.frame(sim.1)

# generate parameters for each simulation. these parameters are constant
# within each simulation but vary across simulations.
sim.1$b0 = 0 + rnorm(nrow(sim.1))
sim.1$b1 = 1 + rnorm(nrow(sim.1))
sim.1$b2 = 3 + rnorm(nrow(sim.1))
sim.1$g1 = 1 + rnorm(nrow(sim.1))
sim.1$g2 = 3 + rnorm(nrow(sim.1))
rownames(sim.1) = paste('tag', sim.1$sim, 'end', sep='.')
sim.0$b0 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'b0']
sim.0$b1 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'b1']
sim.0$b2 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'b2']
sim.0$g1 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'g1']
sim.0$g2 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'g2']

# generate outcome
sim.0$y = with(sim.0, b0 + (b1 * x1) + (b2 * x2) + (g1 * w1) + (g2 * w2) + epsilon)

keeps = c('sim', 'y', 'x1', 'x2', 'w1', 'w2')
sim.2 = sim.0[, keeps]

fold.id.0 = c(t(matrix(c(1, 2), nrow=2, ncol=NOBS/2)))

#===========================================================
# Estimation Start
#===========================================================

#
# Estimation of one simulation at a time
# using linearDML
#
for(i in 1:NSIM){
  if((i %% 10) == 1){
    print(i)
  }

  sim.3 = sim.2[sim.2$sim == sim.1$sim[i], ]

  sim.3$x100 = sim.3$x1 + sim.3$x2
  sim.4 = sim.3[fold.id.0 == 1, ]
  sim.5 = sim.3[fold.id.0 == 2, ]

  # OLS
  model.2 = lm(y ~ x1 + x2 + w1 + w2, data=sim.3)
  coeff.2 = as.data.frame(summary(model.2)$coefficients)
  sim.1[i, c('b1.hat2', 'b2.hat2')] = coeff.2[c('x1', 'x2'), 'Estimate']

  # First Stages
  sim.3[fold.id.0 == 2, 'yt'] = sim.5$y - predict(lm(y ~ w1 + w2, data=sim.4), newdata=sim.5)
  sim.3[fold.id.0 == 2, 'x1t'] = sim.5$x1 - predict(lm(x1 ~ w1 + w2, data=sim.4), newdata=sim.5)
  sim.3[fold.id.0 == 2, 'x2t'] = sim.5$x2 - predict(lm(x2 ~ w1 + w2, data=sim.4), newdata=sim.5)
  sim.3[fold.id.0 == 2, 'x100h'] = predict(lm(x100 ~ w1 + w2, data=sim.4), newdata=sim.5)
  sim.3[fold.id.0 == 1, 'yt'] = sim.4$y - predict(lm(y ~ w1 + w2, data=sim.5), newdata=sim.4)
  sim.3[fold.id.0 == 1, 'x1t'] = sim.4$x1 - predict(lm(x1 ~ w1 + w2, data=sim.5), newdata=sim.4)
  sim.3[fold.id.0 == 1, 'x2t'] = sim.4$x2 - predict(lm(x2 ~ w1 + w2, data=sim.5), newdata=sim.4)
  sim.3[fold.id.0 == 1, 'x100h'] = predict(lm(x100 ~ w1 + w2, data=sim.5), newdata=sim.4)
  sim.3$x100t = sim.3$x100 - sim.3$x100h
  sim.3$x1tt = sim.3$x100t * sim.3$x1
  sim.3$x2tt = sim.3$x100t * sim.3$x2

  # MR
  model.3 = lm(yt ~ x1t + x2t -1, data=sim.3)
  coeff.3 = as.data.frame(summary(model.3)$coefficients)
  sim.1[i, c('b1.hat3', 'b2.hat3')] = coeff.3[c('x1t', 'x2t'), 'Estimate']

  # SR1
  model.4 = lm(yt ~ x1tt + x2tt -1, data=sim.3)
  coeff.4 = as.data.frame(summary(model.4)$coefficients)
  sim.1[i, c('b1.hat4', 'b2.hat4')] = coeff.4[c('x1tt', 'x2tt'), 'Estimate']

  # SR2
  model.5 = lm(yt ~ x1 + x2 + x100h -1, data=sim.3)
  coeff.5 = as.data.frame(summary(model.5)$coefficients)
  sim.1[i, c('b1.hat5', 'b2.hat5')] = coeff.5[c('x1', 'x2'), 'Estimate']


  for(vv in c('mr', 'sr1', 'sr2')){
    model.0 = dml.lm(sim.3, 'y', c('w1', 'w2'), d_vars = c('x1', 'x2'),
                     first_stage_family = 'ols', second_stage_family = vv,
                     foldid=fold.id.0)
    coeff.0 = as.data.frame(summary(model.0$model)$coefficients)

    tmp1 = paste('b1.hat', vv, sep='.')
    tmp2 = paste('b2.hat', vv, sep='.')
    sim.1[i, c(tmp1, tmp2)] = coeff.0[c('x1', 'x2'), 'Estimate']
  }

  # MR - RF prediction
  model.0 = suppressMessages(dml.lm(sim.3, 'y', c('w1', 'w2'), d_vars = c('x1', 'x2'),
                   first_stage_family = 'rf', second_stage_family = 'mr',
                   trees=NTREES, min_n = MIN_N, mtry= 2))
  coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
  sim.1[i, c('b1.hat.rf', 'b2.hat.rf')] = coeff.0[c('x1', 'x2'), 'Estimate']

  # DoubleML
  learner = lrn("regr.ranger", num.trees = NTREES, mtry=2, min.node.size=MIN_N)
  ml_l = learner$clone()
  ml_m = learner$clone()
  obj_dml_data = suppressMessages(DoubleMLData$new(sim.3[, 2:6], y_col="y", d_cols=c('x1', 'x2')))
  dml_plr_obj = suppressMessages(DoubleMLPLR$new(obj_dml_data, ml_l, ml_m, n_folds=2))
  suppressMessages(dml_plr_obj$fit())

  sim.1[i, c('b1.hat.dub', 'b2.hat.dub')] = c(dml_plr_obj$coef)
}

# how does MR differ from long regression?
sim.1$b1.hat.d0 = abs(sim.1$b1.hat2 - sim.1$b1.hat3)
sim.1$b2.hat.d0 = abs(sim.1$b2.hat2 - sim.1$b2.hat3)

# how does MR differ from dml.lm(second_stage_family='mr')
sim.1$b1.hat.dmr = abs(sim.1$b1.hat.mr - sim.1$b1.hat3)
sim.1$b2.hat.dmr = abs(sim.1$b2.hat.mr - sim.1$b2.hat3)

# how does SR1 differ from dml.lm(second_stage_family='sr1')
sim.1$b1.hat.dsr1 = abs(sim.1$b1.hat.sr1 - sim.1$b1.hat4)
sim.1$b2.hat.dsr1 = abs(sim.1$b2.hat.sr1 - sim.1$b2.hat4)

# how does SR2 differ from dml.lm(second_stage_family='sr2')
sim.1$b1.hat.dsr2 = abs(sim.1$b1.hat.sr2 - sim.1$b1.hat5)
sim.1$b2.hat.dsr2 = abs(sim.1$b2.hat.sr2 - sim.1$b2.hat5)

# how does dml.lm(first_stage_family='rf') differ from true parameter?
sim.1$b1.hat.drf = abs(sim.1$b1.hat.rf - sim.1$b1)
sim.1$b2.hat.drf = abs(sim.1$b2.hat.rf - sim.1$b2)

# how does DoubleML with rf differ from true parameter?
sim.1$b1.hat.ddub = abs(sim.1$b1 - sim.1$b1.hat.dub)
sim.1$b2.hat.ddub = abs(sim.1$b2 - sim.1$b2.hat.dub)

# how does dml.lm(first_stage_family='rf') differ from DoubleML?
sim.1$b1.hat.dcmp = abs(sim.1$b1.hat.rf - sim.1$b1.hat.dub)
sim.1$b2.hat.dcmp = abs(sim.1$b2.hat.rf - sim.1$b2.hat.dub)

measure.vars.0 = c('d0', 'dmr', 'dsr1', 'dsr2', 'drf', 'ddub', 'dcmp')
measure.vars.1 = c(paste('b1.hat', measure.vars.0, sep='.'),
                   paste('b2.hat', measure.vars.0, sep='.'))
sim.4 = as.data.frame(melt(sim.1, id.vars='sim',
                           measure.vars=measure.vars.1))
sim.4$one = 1
sim.4$coeff = substr(sim.4$variable, 1, 2)
sim.4$diff.pair = substr(sim.4$variable, 8, 20)

sim.5 = sim.4 %>%
  group_by(diff.pair, coeff) %>%
  summarize(obs=sum(one),
            mean.abs.diff=mean(value),
            rt.mn.sq.diff=sqrt(mean(value^2)),
            .groups='keep')
sim.5 = as.data.frame(sim.5)

print(script.dir)
setwd(script.dir)

write.csv(sim.5, 'test-olsone.csv', row.names=F)


#===========================================================
# Estimation End
#===========================================================

END_TIME = Sys.time()
print(paste('Time elapsed in ', script.name, ':', sep=''))
print(END_TIME - START_TIME)



