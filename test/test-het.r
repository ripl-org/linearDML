START_TIME = Sys.time()

library(dplyr)
library(reshape)
library(stringr)

library(DoubleML)
library(mlr3)
library(mlr3learners)
devtools::load_all()

script.name = basename(sys.frame(1)$ofile)
script.dir = dirname(sys.frame(1)$ofile)

lgr::get_logger("mlr3")$set_threshold("warn")

NTREES = 100
MIN_N = 2
NSIM = 100
NOBS = 1000

sim.0 = expand.grid(id=1:NOBS, sim=1:NSIM)

set.seed(20230216)

sim.0$one = 1
sim.0$c0 = 1

# a0 = 1:99 / 100
# plot(a0, scale(qnorm(a0)))
# plot(a0, scale(qexp(a0)))
# plot(a0, scale(qbeta(a0, 0.1, 0.2)))
# plot(a0, scale(qbinom(a0, size=20, prob=0.08)))
# plot(a0, qunif(a0) + qnorm(a0))
#
# n = 10000
# a1 = rnorm(n) + (rbinom(n, size=2, prob=0.5) * 4)
# hist(a1, breaks=50)
# sd(a1)

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

# delta.0 = rnorm(20)

# primitives - all gaussian
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
sim.0$x1 = as.numeric(sim.0$x1.0 + sim.0$eta1 < 0)

# continuous treatment
sim.0$eta2 = rnorm(nrow(sim.0), mean=sig, sd=sig * 0.9)
sim.0$x2 = sim.0$x2.0 + sim.0$eta2

# # diagnostic
# look.0 = sim.0 %>%
#   group_by(sim, delta1.1, delta2.1, delta1.6, delta2.6,
#            delta1.11, delta2.11, delta1.16, delta2.16) %>%
#   summarize(obs=sum(one),
#             x1.w1=cor(x1, w1),
#             x2.w1=cor(x2, w1),
#             x1.w6=cor(x1, w6),
#             x2.w6=cor(x2, w6),
#             x1.w11=cor(x1, w11),
#             x2.w11=cor(x2, w11),
#             x1.w16=cor(x1, w16),
#             x2.w16=cor(x2, w16),
#             .groups='keep')
#
# # diagnostic
# look.1 = sim.0[sim.0$sim == 5, ]
# # plot(look.1$w6, look.1$u6)
# hist(look.1$y1)


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
sim.1$beta1.10 = rnorm(nrow(sim.1))
sim.0$beta1.10 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'beta1.10']

sim.1$beta2.10 = rnorm(nrow(sim.1))
sim.0$beta2.10 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'beta2.10']

sim.1$beta1.20 = rnorm(nrow(sim.1))
sim.0$beta1.20 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'beta1.20']

sim.1$beta2.20 = rnorm(nrow(sim.1))
sim.0$beta2.20 = sim.1[paste('tag', sim.0$sim, 'end', sep='.'), 'beta2.20']

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

# 5. no heterogeneity, x1 and x2
# 6. yes heterogeneity, x1 and x2


keeps = c('sim', 'c0', paste('y', 1:4, sep=''), 'x1', 'x2', paste('w', 1:20, sep=''))
sim.2 = sim.0[, keeps]

# fold.id.0 = c(t(matrix(c(1, 2), nrow=2, ncol=NOBS/2)))

# stop('here')

#===========================================================
# Estimation Start
#===========================================================


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
w.all = paste('w', 1:20, sep='')

sim.2$x1_w10 = with(sim.2, x1 * w10)
sim.2$x1_w20 = with(sim.2, x1 * w20)
sim.2$x2_w10 = with(sim.2, x2 * w10)
sim.2$x2_w20 = with(sim.2, x2 * w20)


#
# Estimation of one simulation at a time
# using linearDML
#
for(i in 1:NSIM){
  if((i %% 10) == 1){
    print(i)
  }

  sim.3 = sim.2[sim.2$sim == sim.1$sim[i], ]

  for(j in 1:4){
    h.vars = get(paste('h.vars', j, sep='.'))
    y.curr = paste('y', j, sep='')
    x.curr = paste('x', ceiling(j / 2), sep='')
    x.curr.10 = paste(x.curr, 'w10', sep='_')
    x.curr.20 = paste(x.curr, 'w20', sep='_')
    d.vars = h.vars$fxd.name
    if((j %% 2) == 0){
      jj = 1:3
    }else{
      jj = 1
    }
    ajj = c(x.curr, x.curr.10, x.curr.20)[jj]

    #
    # DML with OLS prediction
    #
    TIME_00 = Sys.time()
    model.0 = suppressMessages(dml.lm(sim.3, y.curr, w.all, h_vars=h.vars,
                                      first_stage_family = 'ols', second_stage_family = 'mr'))
    coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
    colnames(coeff.0) = c('b', 'se', 't', 'p')
    aj = paste(c('a.hat.ols', 'b10.hat.ols', 'b20.hat.ols'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'b']
    aj = paste(c('a.se.ols', 'b10.se.ols', 'b20.se.ols'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'se']
    TIME_01 = Sys.time()
    sim.1[i, paste('time.ols', j, sep='.')] = TIME_01 - TIME_00

    #
    # DML with RF prediction
    #
    TIME_00 = Sys.time()
    model.0 = suppressMessages(dml.lm(sim.3, y.curr, w.all, h_vars=h.vars,
                                      first_stage_family = 'rf', second_stage_family = 'mr',
                                      trees = NTREES, mtry=10, min_n=MIN_N))
    coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
    colnames(coeff.0) = c('b', 'se', 't', 'p')
    aj = paste(c('a.hat.rf', 'b10.hat.rf', 'b20.hat.rf'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'b']
    aj = paste(c('a.se.rf', 'b10.se.rf', 'b20.se.rf'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'se']
    TIME_01 = Sys.time()
    sim.1[i, paste('time.rf', j, sep='.')] = TIME_01 - TIME_00

    #
    # DML with RF prediction, not calculated heterogeneously
    #
    TIME_00 = Sys.time()
    model.0 = suppressMessages(dml.lm(sim.3, y.curr, w.all, d_vars=d.vars,
                                      first_stage_family = 'rf', second_stage_family = 'mr',
                                      trees = NTREES, mtry=10, min_n=MIN_N))
    coeff.0 = as.data.frame(summary(model.0$model)$coefficients)
    colnames(coeff.0) = c('b', 'se', 't', 'p')
    aj = paste(c('a.hat.rf2', 'b10.hat.rf2', 'b20.hat.rf2'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'b']
    aj = paste(c('a.se.rf2', 'b10.se.rf2', 'b20.se.rf2'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'se']
    TIME_01 = Sys.time()
    sim.1[i, paste('time.rf2', j, sep='.')] = TIME_01 - TIME_00



    #
    # DML from DoubleML
    # 2023-02-16 AYSM: oh actually maybe this doesn't allow us to see
    # heterogeneous treatment effects at all, it just relaxes the assumption
    # so that there *could plausibly be* heterogeneous treatment effects
    #
    # 2023-02-21 AYSM: we can always just implement it in the partially linear
    # model including the heterogeneous treatments as treatments, but will be
    # different

    TIME_00 = Sys.time()
    learner = lrn("regr.ranger", num.trees = NTREES, mtry=10, min.node.size=MIN_N)
    ml_l = learner$clone()
    ml_m = learner$clone()
    obj_dml_data = suppressMessages(DoubleMLData$new(sim.3[, c(y.curr, d.vars, w.all)],
                                                     y_col=y.curr, d_cols=d.vars))
    dml_plr_obj = suppressMessages(DoubleMLPLR$new(obj_dml_data, ml_l, ml_m, n_folds=2))
    suppressMessages(dml_plr_obj$fit())

    # ml_m_2 = lrn("classif.ranger", num.trees = NTREES, mtry=10, min.node.size=MIN_N)
    # dml_irm_obj = DoubleMLIRM$new(obj_dml_data, ml_l, ml_m_2, n_folds=2)
    # dml_irm_obj$fit(store_predictions = T)

    # dml_irm_obj$coef
    # dml_irm_obj$print()
    #
    # dml_irm_obj$summary()
    # dml_irm_obj$fit()
    #
    # sim.3$g0 = dml_irm_obj$predictions$ml_g0
    # sim.3$g1 = dml_irm_obj$predictions$ml_g1
    # sim.3$te = with(sim.3, g1 - g0)
    # plot(sim.3$w10, sim.3$te)
    #
    # mean(sim.3[sim.3$x1 >= 0, 'te'])


    invisible(capture.output(coeff.0 <- as.data.frame(as.matrix(dml_plr_obj$summary()))))
    colnames(coeff.0) = c('b', 'se', 't', 'p')
    aj = paste(c('a.hat.dub', 'b10.hat.dub', 'b20.hat.dub'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'b']
    aj = paste(c('a.se.dub', 'b10.se.dub', 'b20.se.dub'), j, sep='.')[jj]
    sim.1[i, aj] = coeff.0[ajj, 'se']
    TIME_01 = Sys.time()
    sim.1[i, paste('time.dub', j, sep='.')] = TIME_01 - TIME_00

  }
}



measure.vars.0 = c()
for(vv in c('ols', 'rf', 'rf2', 'dub')){
  for(j in 1:4){
    sim.1[, paste('da', vv, j, sep='.')] =
      abs(sim.1$alpha1 - sim.1[, paste('a.hat', vv, j, sep='.')])

    measure.vars.0 = c(measure.vars.0, paste('da', vv, j, sep='.'))
  }

  sim.1[, paste('db10', vv, 2, sep='.')] = abs(sim.1$beta1.10 - sim.1[, paste('b10.hat', vv, 2, sep='.')])
  sim.1[, paste('db20', vv, 2, sep='.')] = abs(sim.1$beta1.20 - sim.1[, paste('b20.hat', vv, 2, sep='.')])
  sim.1[, paste('db10', vv, 4, sep='.')] = abs(sim.1$beta2.10 - sim.1[, paste('b10.hat', vv, 4, sep='.')])
  sim.1[, paste('db20', vv, 4, sep='.')] = abs(sim.1$beta2.20 - sim.1[, paste('b20.hat', vv, 4, sep='.')])

  measure.vars.0 = c(measure.vars.0, paste(c('db10', 'db20'), vv, 2, sep='.'),
                     paste(c('db10', 'db20'), vv, 4, sep='.'))
}



sim.4 = as.data.frame(melt(sim.1, id.vars='sim',
                           measure.vars=measure.vars.0))

sim.4$one = 1
tmp1 = str_match(sim.4$variable, '(.*)\\.(.*)\\.(.*)')
sim.4$coeff = substr(tmp1[, 2], 2, 10)
sim.4$s1 = tmp1[, 3]
sim.4$y.cat = tmp1[, 4]

sim.5 = sim.4 %>%
  group_by(y.cat, coeff, s1, var.name=variable) %>%
  summarize(obs=sum(one),
            mean.abs.diff=mean(value),
            rt.mn.sq.diff=sqrt(mean(value^2)),
            .groups='keep')
sim.5 = as.data.frame(sim.5)


measure.vars.0 = colnames(sim.1)[grepl('se', colnames(sim.1))]
sim.6 = as.data.frame(melt(sim.1, id.vars='sim',
                            measure.vars=measure.vars.0))

# a.se.dub.j
tmp1 = str_match(sim.6$variable, '(.*)\\.(.*)\\.(.*)\\.(.*)')
sim.6$one = 1
sim.6$coeff = tmp1[, 2]
# tmp1 = nchar(as.character(sim.6$variable))
sim.6$s1 = tmp1[, 4]
sim.6$y.cat = tmp1[, 5]

sim.7 = sim.6 %>%
  group_by(y.cat, coeff, s1) %>%
  summarize(se.mean=mean(value),
            .groups='keep')
sim.7 = as.data.frame(sim.7)

sim.8 = merge(sim.5, sim.7, by=c('y.cat', 'coeff', 's1'))


measure.vars.0 = colnames(sim.1)[grepl('time', colnames(sim.1))]
time.0 = as.data.frame(melt(sim.1, id.vars='sim',
                           measure.vars=measure.vars.0))
time.0$one = 1
tmp1 = nchar(as.character(time.0$variable))
time.0$s1 = substr(time.0$variable, 6, tmp1 - 2)
time.0$y.cat = substr(time.0$variable, tmp1, tmp1)

time.1 = time.0 %>%
  group_by(y.cat, s1, var.name=variable) %>%
  summarize(obs=sum(one),
            mean.time=mean(value),
            .groups='keep')
time.1 = as.data.frame(time.1)



print(script.dir)
setwd(script.dir)

write.csv(sim.8, 'test-hetcoef.csv', row.names=F)
write.csv(time.1, 'test-hettime.csv', row.names=F)


#===========================================================
# Estimation End
#===========================================================

END_TIME = Sys.time()
print(paste('Time elapsed in ', script.name, ':', sep=''))
print(END_TIME - START_TIME)



