START_TIME = Sys.time()

library(dplyr)
library(reshape)

devtools::load_all()

script.name = basename(sys.frame(1)$ofile)


NSIM = 200
NOBS = 500
sim.0 = expand.grid(id=1:NOBS, sim=1:NSIM)

set.seed(20220822)

sim.0$one = 1

# RHO = 1

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


W = as.matrix(sim.0[, c('c0', 'w1', 'w2')])
X = as.matrix(sim.0[, c('x1', 'x2')])

# as stated above, the covariance matrix between the x's and
# the w's is full rank.
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

# generate parameters for each simulation
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

# split variables by simulation
il = sim.1$ind.lo
ih = sim.1$ind.hi
Y_s = lapply(seq_len(NSIM), function(x) as.matrix(
  sim.0[il[x]:ih[x], c('y')]))
C_s = lapply(seq_len(NSIM), function(x) as.matrix(
  sim.0[il[x]:ih[x], c('c0')]))
D_s = lapply(seq_len(NSIM), function(x) as.matrix(
  sim.0[il[x]:ih[x], c('x1', 'x2')]))
Db_s = lapply(seq_len(NSIM), function(x) as.matrix(
  sim.0[il[x]:ih[x], c('xbar')]))
W_s = lapply(seq_len(NSIM), function(x) as.matrix(
  sim.0[il[x]:ih[x], c('w1', 'w2')]))

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

  for(vv in c('mr', 'sr1', 'sr2')){
    model.0 = dml.lm(sim.3, 'y', c('w1', 'w2'), d_vars = c('x1', 'x2'),
                     first_stage_family = 'ols', second_stage_family = vv,
                     foldid=fold.id.0)
    coeff.0 = as.data.frame(summary(model.0)$coefficients)

    tmp1 = paste('b1.hat', vv, sep='.')
    tmp2 = paste('b2.hat', vv, sep='.')
    sim.1[i, c(tmp1, tmp2)] = coeff.0[c('x1', 'x2'), 'Estimate']
  }

  # stop('here 126')
}


#
# Mass Estimation
#

# long regression
X_s = Map(function(cc, dd, ww) cbind(cc, dd, ww), cc=C_s, dd=D_s, ww=W_s)
B_1 = Map(function(x, y) c(solve(t(x) %*% x) %*% (t(x) %*% y)), x=X_s, y=Y_s)
B_all = do.call(rbind, B_1)
sim.1[, c('b0.hat2', 'b1.hat2', 'b2.hat2', 'g1.hat2', 'g2.hat2')] = B_all


cross.resid <- function(y, x, foldid){
  ya = y[foldid == 1, 1:ncol(y), drop=F]
  yb = y[foldid == 2, 1:ncol(y), drop=F]
  xa = x[foldid == 1, 1:ncol(x), drop=F]
  xb = x[foldid == 2, 1:ncol(x), drop=F]

  ba = solve(t(xa) %*% xa) %*% (t(xa) %*% ya)
  bb = solve(t(xb) %*% xb) %*% (t(xb) %*% yb)

  ret = y
  ret[foldid == 1, 1:ncol(y)] = ya - (xa %*% bb)
  ret[foldid == 2, 1:ncol(y)] = yb - (xb %*% ba)
  return(ret)
}


# dml
X_s = D_s
Z_s = Map(function(cc, ww) cbind(cc, ww), cc=C_s, ww=W_s)
Yt_s = Map(function(w, y) (cross.resid(y, w, fold.id.0)),
           w=Z_s, y=Y_s)
Xt_s = Map(function(w, x) (cross.resid(x, w, fold.id.0)),
           w=Z_s, x=X_s)
B_2 = Map(function(x, y) c(solve(t(x) %*% x) %*% (t(x) %*% y)), x=Xt_s, y=Yt_s)
B_all = do.call(rbind, B_2)
sim.1[, c('b1.hat3', 'b2.hat3')] = B_all



# pseudo het dml
Z_s = Map(function(cc, ww) cbind(cc, ww), cc=C_s, ww=W_s)
Yt_s = Map(function(w, y) (cross.resid(y, w, fold.id.0)),
           w=Z_s, y=Y_s)
Dbt_s = Map(function(w, d) (cross.resid(d, w, fold.id.0)),
            w=Z_s, d=Db_s)
X_s = Map(function(d, dbt, cc) cbind((d * cbind(dbt, dbt))),
          d=D_s, dbt=Dbt_s, cc=C_s)
B_3 = Map(function(x, y) c(solve(t(x) %*% x) %*% (t(x) %*% y)), x=X_s, y=Yt_s)
B_all = do.call(rbind, B_3)
sim.1[, c('b1.hat4', 'b2.hat4')] = B_all



# pseudo het dml - fixed
Dbt_s = Map(function(w, d) (d - cross.resid(d, w, fold.id.0)),
            w=Z_s, d=Db_s)
X_s = Map(function(d, dbt, cc) cbind(dbt, d),
          d=D_s, dbt=Dbt_s, cc=C_s)
B_3 = Map(function(x, y) c(solve(t(x) %*% x) %*% (t(x) %*% y)), x=X_s, y=Yt_s)
B_all = do.call(rbind, B_3)
sim.1[, c('b100.hat5', 'b1.hat5', 'b2.hat5')] = B_all


sim.1$b1.hat.d0 = abs(sim.1$b1.hat2 - sim.1$b1.hat3)
sim.1$b2.hat.d0 = abs(sim.1$b2.hat2 - sim.1$b2.hat3)

sim.1$b1.hat.dmr = abs(sim.1$b1.hat.mr - sim.1$b1.hat3)
sim.1$b2.hat.dmr = abs(sim.1$b2.hat.mr - sim.1$b2.hat3)

sim.1$b1.hat.dsr1 = abs(sim.1$b1.hat.sr1 - sim.1$b1.hat4)
sim.1$b2.hat.dsr1 = abs(sim.1$b2.hat.sr1 - sim.1$b2.hat4)

sim.1$b1.hat.dsr2 = abs(sim.1$b1.hat.sr2 - sim.1$b1.hat5)
sim.1$b2.hat.dsr2 = abs(sim.1$b2.hat.sr2 - sim.1$b2.hat5)

sim.1$one = 1
sim.4 = sim.1 %>%
  summarize(nsim=sum(one),
            b1.hat.d0=sum(b1.hat.d0),
            b2.hat.d0=sum(b2.hat.d0),
            b1.hat.dmr=sum(b1.hat.dmr),
            b2.hat.dmr=sum(b2.hat.dmr),
            b1.hat.dsr1=sum(b1.hat.dsr1),
            b2.hat.dsr1=sum(b2.hat.dsr1),
            b1.hat.dsr2=sum(b1.hat.dsr2),
            b2.hat.dsr2=sum(b2.hat.dsr2))



#===========================================================
# Estimation End
#===========================================================

END_TIME = Sys.time()
print(paste('Time elapsed in ', script.name, ':', sep=''))
print(END_TIME - START_TIME)



