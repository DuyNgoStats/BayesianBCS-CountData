library(credsubs)
library(rjags)
library(mpath)
library(pscl)
library(VGAM)

source('BCS.R')
source('sim_jags.R')


##########################################################################
##########################################################################
##########################################################################
n=250     ###70, 150, 250
cred.level=0.8
nsim=1000
nIter=10000
zeta.diff = 0
cs.spikelab = matrix(0,nsim,5)
colnames(cs.spikelab) = c('Total Coverage', 'Pair Size', 'Sensitivity', "Specificity", 'MSE')
cs.pointwise = cs.hh = cs.spikelab

results= matrix(0,nsim,13)

beta = c(0.7, -0.1, 0.1,  0.4, 0, 0.01)

alpha = c(-1, 0, -0.2, -0.8, 0, 0.1) #logit 

#alpha = c(-1.1, 0.1, 0, -0.9, 0.01, 0) #cloglog

#alpha = c(-0.9, -0.2, 0, 0.2, 0.1, 0) #probit


for(i in 1:nsim){
  print(paste('Iter ', i))
  
  X = matrix(0, n, 6)
  X[,1] = 1 ## intercept
  X[,2] = rbinom(n, 1, 0.6) 
  X[,3] = rnorm(n, 0, 1) 
  X[,4] = rbinom(n,1,0.5) ##treatment
  X[,5] = X[,2] * X[, 4]
  X[,6] = X[,3] * X[, 4]
  
  Z = matrix(0, n, 6)
  Z[,1] = 1 ## intercept
  Z[,2] = rbinom(n, 1, 0.3)  
  Z[,3] = X[,3]
  Z[,4] = X[,4] # treatment
  Z[,5] = Z[,2] * Z[,4]
  Z[,6] = Z[,3] * Z[,4]
  
  data=simulateData(n, X, Z, beta, alpha)
  
  datajags = list(y = data$y, X = X, Z=Z)
  
  JagResults=PoissonJag(textConnection(zipModel_spikeslab),datajags,nIter)
  JagResults2=PoissonJag(textConnection(zipModel_horseshoe),datajags,nIter)
 
  ############################################################
  ############################################################
   ### design matrix
  # Construct grid on predictive covariate region of interest
  designTrt = expand.grid(x1 = range(X[,2]),
                          x2 = seq(min(X[,3]), max(X[,3]), length.out=10),
                          x3 = c(1),
                          z1 = range(Z[,2]),
                          z2 = seq(min(Z[,3]), max(Z[,3]), length.out=10))
                          
  
  designCtrl = expand.grid(x1 = range(X[,2]),
                           x2 = seq(min(X[,3]), max(X[,3]), length.out=10),
                           x3 =  c(0),
                           z1 = range(Z[,2]),
                           z2 = seq(min(Z[,3]), max(Z[,3]), length.out=10))
  
  
  pte.spikeslab = ComputePTE(designTrt, designCtrl, coef.zero=JagResults$pos.logit, 
             coef.pois=JagResults$pos.log, alpha, beta)
  
  pte.hh = ComputePTE(designTrt, designCtrl, coef.zero=JagResults2$pos.logit, 
                             coef.pois=JagResults2$pos.log, alpha, beta)
  
  ### Compute estimated personalized treatment effect
  real.safety.diff  = pte.spikeslab$real.safety.diff
  
  #### Compute the true personalized treatment effect
  real.subgroup.fromDesign.diff <- real.safety.diff  > zeta.diff
  
  cs.spikelab[i,] = CSSafety(est.safety= pte.spikeslab$safety.diff, real.safety = real.safety.diff, 
                         real.subgroup.fromDesign=real.subgroup.fromDesign.diff, 
                         zeta=zeta.diff, cred.level=cred.level)
  cs.pointwise[i,] = CSSafetyPointwise(safety.diff = pte.spikeslab$safety.diff, df = nrow(pte.spikeslab$safety.diff)-1, 
                                       real.safety.diff, real.subgroup.fromDesign.diff,
                                       zeta=zeta.diff, cred.level=cred.level)
  cs.hh[i,] = CSSafety(est.safety= pte.hh$safety.diff, real.safety = real.safety.diff, 
                       real.subgroup.fromDesign=real.subgroup.fromDesign.diff, 
                       zeta=zeta.diff, cred.level=cred.level)
  #results[i,]=JagResults$posteriorSummary[1:13]
}



colMeans(cs.spikelab)
colMeans(cs.pointwise)
colMeans(cs.hh)

