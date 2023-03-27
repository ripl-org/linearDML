START_TIME = Sys.time()

library(dplyr)
library(reshape)
library(stringr)

# Description:
# This file tests the bias of the dml.lm() for each second_stage_family
# option against the theoretically predicted bias.

# Load the package
devtools::load_all()
# Find the current script
script.name = basename(sys.frame(1)$ofile)
script.dir = dirname(sys.frame(1)$ofile)

# Simulation parameters
NSIM = 50
NOBS = 10000

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

b0 = 0 + rnorm(1)
b1 = 1 + rnorm(1)
b2 = 3 + rnorm(1)
g1 = 1 + rnorm(1)
g2 = 3 + rnorm(1)

# generate parameters for each simulation. these parameters are constant
# within each simulation and across simulations.
sim.1$b0 = b0
sim.1$b1 = b1
sim.1$b2 = b2
sim.1$g1 = g1
sim.1$g2 = g2
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
    print(paste(Sys.time(), i))
  }

  sim.3 = sim.2[sim.2$sim == sim.1$sim[i], ]

  for(vv in c('mr', 'sr1', 'sr2')){
    for(mm in c('ols', 'rf')){
      model.0 = dml.lm(sim.3, 'y', c('w1', 'w2'), d_vars = c('x1', 'x2'),
                       first_stage_family = mm, second_stage_family = vv,
                       foldid=fold.id.0)
      coeff.0 = as.data.frame(summary(model.0)$coefficients)

      tmp1 = paste('b1.hat', vv, mm, sep='.')
      tmp2 = paste('b2.hat', vv, mm, sep='.')
      sim.1[i, c(tmp1, tmp2)] = coeff.0[c('x1', 'x2'), 'Estimate']
    }
  }
}


measure.vars.0 = names(sim.1)[grepl('b[12].hat', names(sim.1))]
sim.11 = as.data.frame(melt(sim.1, id.vars=c('sim'),
                            measure.vars=measure.vars.0))
names(sim.11) = c('sim', 'variable', 'estimate')
sim.11$var.name = substr(sim.11$variable, 1, 2)
sim.11$method = substr(sim.11$variable, 8, 20)
sim.11[, c('s2', 's1')] = str_split_fixed(sim.11$method, '\\.', n=2)
sim.11$true.value = with(sim.11, ifelse(var.name == 'b1', b1, b2))
sim.11$diff = with(sim.11, estimate - true.value)

# Get bias among simulations as average difference with true value
sim.11$one = 1
sim.5 = sim.11 %>%
  group_by(variable, var.name, method, s1, s2) %>%
  summarize(nsim=sum(one),
            dmn=mean(diff),
            dsd=sd(diff),
            .groups='keep')
sim.5 = as.data.frame(sim.5)
sim.6 = as.data.frame(melt(sim.5, id.vars=c('var.name', 'method', 's1', 's2'),
                           measure.vars=c('dmn', 'dsd')))
sim.7 = as.data.frame(cast(sim.6, s2 + var.name ~ s1 + variable, value='dmn'))

#===========================================================
# Estimation End
#===========================================================

#
# Now, calculate the theoretical bias of SR1 and SR2
#

PSI_0 = as.matrix(sim.0[, c('c0', 'w1', 'w2')]) %*% Delta_0

# SR2
Z_0 = PSI_0 %*% rbind(1, 1)
D_0 = as.matrix(sim.0[, c('x1', 'x2')])
THETA_0 = solve(t(Z_0) %*% Z_0) %*% (t(Z_0) %*% D_0)
Z_THETA = Z_0 %*% THETA_0

ZT2 = t(Z_THETA) %*% Z_THETA / (NSIM * NOBS)
PSI2 = t(PSI_0) %*% PSI_0 / (NSIM * NOBS)
D2 = t(D_0) %*% D_0 / (NSIM * NOBS)

f1 = (D2) - (ZT2)
f2 = (PSI2) - (ZT2)
f3 = solve(f1) %*% f2
f4 = -f3 %*% rbind(b1, b2)

sim.7[5:6, 'approx.bias'] = f4[1:2, 1]

# SR1
f5 = D_0 %*% rbind(b1, b2)
f6 = Z_0 * f5
f7 = PSI_0 %*% rbind(b1, b2)
f8 = f6 - f7
D_TIL = D_0 - cbind(Z_0, Z_0) * D_0
f9 = t(D_TIL) %*% f8 / (NSIM * NOBS)
f10 = t(D_TIL) %*% D_TIL / (NSIM * NOBS)
f11 = solve(f10) %*% f9

sim.7[3:4, 'approx.bias'] = f11[1:2, 1]


print(script.dir)
setwd(script.dir)

write.csv(sim.7, 'test-bias.csv', row.names=F)

END_TIME = Sys.time()
print(paste('Time elapsed in ', script.name, ':', sep=''))
print(END_TIME - START_TIME)



