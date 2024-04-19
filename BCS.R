library(credsubs)
library(rjags)
library(mpath)
library(pscl)
library(VGAM)

### n: number of subjects.
### X, beta: covariate and coef of Poisson component.
### Z, alpha: covariate and coef of zero component.
simulateData = function(n, X, Z, beta, alpha){
  
  #Mean of the Poisson component.
  mu = exp(X%*%beta)
  
  #Product of Z matrix and coefficients for zero component.
  #theta = probability of rv from zero component of rv.
  mu.z = Z%*%alpha
  
  theta = logitlink(mu.z, inverse = TRUE)
  #theta = clogloglink(mu.z, inverse = TRUE)
  #theta = probitlink(mu.z, inverse = TRUE)
  
  
  zero <- rbinom(n, size=1, prob=theta)
  #indicator variable (if z=0 means poison distribution,z=1 then zero component)
  
  y <- rpois(n, lambda=mu)*(1-zero)
  
  return(list(y=y, X=X, Z=Z, beta=beta, alpha=alpha))
}

################################################################
#######################################################################
PoissonJag = function(jagmodel, data, nIter){
  #Z matrix is the design matrix for zero component of ZIP
  #X matrix is the design matrix for poisson  component of ZIP
  y = data$y
  n = length(y)
  Z = data$Z
  X = data$X
  nKx=ncol(X)
  nKz=ncol(Z)
  
  params = c('alpha','beta')
  
  datalist1 <- list(X=X, Z=Z,y=y,n=n,nKx=nKx,nKz=nKz)
  
  ###JAGS for zipmodel1
  jagsmodel1 <- jags.model(jagmodel, data=datalist1, quiet=TRUE)
  
  
  ###Posterior Sample for zipmodel1
  posteriorSample<- coda.samples(jagsmodel1, params , n.thin=20, 
                                 n.chains=3, n.burnin=500, 
                                 n.iter=nIter, progress.bar='none')
  posteriorSample=posteriorSample[[1]]
  
  return(list(pos.logit = as.matrix(posteriorSample[,1:nKz]),
              pos.log = as.matrix(posteriorSample[,(nKz+1):(nKz + nKx)]),
              posteriorSummary= round(summary(posteriorSample)[[1]][,1],4))
  )
  
}

################################################################
################################################################
CalRegression = function(design, coef.zero, coef.pois, real.coef.zero, real.coef.pois){
  nIter = nrow(coef.zero)# number of iterations
  nsubj = nrow(design)# different characteristics of patient
  #design=designCtrl
  X.design = cbind(intercept = rep(1, nsubj),
                   x1=design$x1,
                   x2=design$x2,
                   x3=design$x3,
                   x4=design$x1 * design$x3,
                   x5=design$x2 * design$x3
  )
  
  Z.design = cbind(intercept = rep(1, nsubj),
                   z1 = design$z1,
                   z2 = design$z2,
                   z3 = design$x3,
                   z4 = design$z1 * design$x3,
                   z5 = design$z2 * design$x3
  )
  
  ##empty matrices
  e.zero = matrix(0, nIter, nsubj)#(1-theta) = probability rv is poisson
  e.pois = matrix(0, nIter, nsubj)#(mu) mean of the poisson.
  
  for(isubj in 1:nsubj){
    e.zero[,isubj] = as.vector(1/(1+exp(Z.design[isubj,]%*%t(coef.zero))))
    e.pois[,isubj] = as.vector(exp(X.design[isubj,]%*%t(coef.pois)))
  }
  
  ###Compute the truth
  real.zero = as.vector(1/(1+exp(Z.design%*%real.coef.zero)))#(1-theta) = true probability rv is poisson
  real.pois = as.vector(exp(X.design%*%real.coef.pois))#(mu) true mean of the poisson.
  
  
  return(list(e.zero=e.zero, e.pois=e.pois, real.zero = real.zero, real.pois = real.pois))
}

################################################################
################################################################
ComputePTE = function(designTrt, designCtrl, coef.zero, coef.pois, alpha, beta){
  ####Expected mean of response for treatment group
  e.trt = CalRegression(design=designTrt, coef.zero = coef.zero, 
                        coef.pois=coef.pois, real.coef.zero = alpha,
                        real.coef.pois = beta)
  
  ####Expected mean of response for treatment group
  e.ctrl = CalRegression(design=designCtrl, coef.zero = coef.zero, 
                         coef.pois=coef.pois, real.coef.zero = alpha,
                         real.coef.pois = beta)
  
  nsubj = nrow(designTrt)
  nIter = nrow(coef.pois)
  safety.diff = matrix(0, nIter, nsubj)
  
  ##############################################
  nsubj = nrow(designTrt )
  
  for(isubj in 1:nsubj){
    safety.diff[,isubj] = (e.trt$e.zero[,isubj])*e.trt$e.pois[,isubj] - 
      (e.ctrl$e.zero[,isubj])*e.ctrl$e.pois[,isubj]
  }
  
  real.safety.diff  = (e.trt$real.zero)*e.trt$real.pois - (e.ctrl$real.zero)*e.ctrl$real.pois
  
  return(list(safety.diff=safety.diff, real.safety.diff=real.safety.diff))
}

################################################################
################################################################
CSSafetyPointwise <- function(safety.diff, df, real.safety.diff,
                              real.subgroup.fromDesign,
                              zeta= 1, cred.level=0.8) {
  
  se = apply(safety.diff, 2, sd)
  mdq = qt(cred.level, nrow(safety.diff) - 1)
  post.mean = apply(safety.diff, 2, mean)
  
  exclusive <- as.vector(post.mean - mdq * se > zeta)
  inclusive <- as.vector(post.mean + mdq * se >= zeta)
  
  cs = list(exclusive = exclusive, inclusive = inclusive)
  
  # Credible subgroup diagnostics
  
  mse.effect <- mean((apply(safety.diff, 2, mean) - real.safety.diff) ^ 2)
  ##D: exclusive, S: inclusive
  csp.size     <- pair.size(cs$exclusive, cs$inclusive)
  summ.stats.D <- diagnostic.test (cs$exclusive, real.subgroup.fromDesign)
  summ.stats.S <- diagnostic.test (cs$inclusive, real.subgroup.fromDesign)
  
  success.D <- all(real.subgroup.fromDesign[cs$exclusive]) # D is exclusive
  success.S <- all(cs$inclusive[real.subgroup.fromDesign]) # S is inclusive
  success <- success.D && success.S
  
  
  return(c(success, csp.size, summ.stats.D, mse.effect)) 
}

################################################################
################################################################
CSSafety = function(est.safety, real.safety, real.subgroup.fromDesign, zeta=1, cred.level=0.8){
  cs = credsubs(params = est.safety, design = NULL, threshold= zeta, cred.level = cred.level)
  
  # Credible subgroup diagnostics
  ## MSE
  mse.effect <- mean((apply(est.safety, 2, mean) - real.safety) ^ 2)
  ##D: exclusive, S: inclusive
  csp.size     <- pair.size(cs$exclusive, cs$inclusive)
  summ.stats.D <- diagnostic.test (cs$exclusive, real.subgroup.fromDesign)
  summ.stats.S <- diagnostic.test (cs$inclusive, real.subgroup.fromDesign)
  #avg.over     <- average.overestimate(cs$exclusive, real.safety, zeta)
  #max.over     <- max(maximum.overestimate(cs$exclusive, real.safety, zeta), 0)
  # Check credible subgroup coverage
  success.D <- all(real.subgroup.fromDesign[cs$exclusive]) # D is exclusive
  success.S <- all(cs$inclusive[real.subgroup.fromDesign]) # S is inclusive
  success <- success.D && success.S
  return(c(success, csp.size, summ.stats.D, mse.effect)) 
}


################################################################
################################################################

pair.size <- function(D,S) {
 
  return(mean(S & !D))
}

################################################################
################################################################

diagnostic.test <- function(est.subset, true.subset) {
 
  true.positive <- true.subset
  true.negative <- !true.positive
  test.positive <- est.subset
  test.negative <- !est.subset
  
sensitivity <- sum(test.positive & true.positive) / sum(true.positive)
#sensitivity is the proportion of true positive
specificity <- sum(test.negative & true.negative) / sum(true.negative)
#specificity is the proportion of true negative
  
  return(c(sensitivity=sensitivity,specificity=specificity))
}


