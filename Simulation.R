source('BCS.R')
source('sim_jags.R')

##############################################################################
##############################################################################
###  TRUE CATE 
beta = matrix(c(1, 0, -0.3,  0, 0, 0,
                1, 0, -0.3,  0.1, 0, 0.5,
                1, 0, -0.3,  0.5, 0, 0.5,
                1, 0, -0.3,  0.7, 0, 0.5), 4, 6, byrow = TRUE)
alpha = matrix(c(-1, 0, 0.5, 0, 0, 0,
                 -1, 0, 0.5, -0.5, 0, -0.1), 2, 6, byrow = TRUE)


pte.null = ComputeTruePTE(beta[1,], alpha[1,], seq(-4, 4, length=100))
pte.small = ComputeTruePTE(beta[2, ], alpha[2,], seq(-4, 4, length=100))
pte.moderate = ComputeTruePTE(beta[3,], alpha[2,], seq(-4, 4, length=100))
pte.large = ComputeTruePTE(beta[4,], alpha[2,], seq(-4, 4, length=100))


plot(seq(-4, 4, length=100), pte.null, ylab = 'CATE', xlab = expression(x[2]), 
     type = 'l', 
     ylim = c(min(c(pte.null, pte.small, pte.moderate, pte.large)) - 0.5, max(c(pte.null, pte.small, pte.moderate, pte.large)) + 0.5), lty = 1,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd=2)
lines(seq(-4, 4, length=100), pte.small, lty = 2, lwd=2)
lines(seq(-4, 4, length=100), pte.moderate, lty = 3, lwd=2)
lines(seq(-4, 4, length=100), pte.large, lty = 4, lwd=2)
legend(0.5, -0.5, legend=c( '(S1) Null', 
                            '(S2) Small effect',
                            '(S3) Moderate effect',
                            '(S4) Large effect'), lwd=3,
       lty=1:4, cex=1.55,title = '', text.font=2, box.lty=0)


#############################################################################
#############################################################################
### Null case
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0,
            x2 = c(-4, 4), link = 'logit', beta = beta[1,], alpha = alpha[1,])


### Small Effect
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0, 
            x2 = c(-4, 4), link = 'logit', beta = beta[2,], alpha = alpha[2,])


### Moderate effect
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0, 
            x2 = c(-4, 4), link = 'logit', beta = beta[3,], alpha = alpha[2,])


### Large effect
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0, 
            x2 = c(-4, 4), link = 'logit', beta = beta[4,], alpha = alpha[2,])


### Full case
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0, 
            x2 = c(0, 4), link = 'logit', beta = beta[3,], alpha = alpha[2,])


########################################################################
########################################################################
### Different link functions
### Moderate effect
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0, 
            x2 = c(-4, 4), link = 'cloglog', beta = beta[3,], alpha = alpha[2,])


### Moderate effect
RunSimulate(sample.size = c(50, 100, 500), nsim = 1000, nIter = 10000, 
            cred.level = 0.8, zeta.diff = 0, 
            x2 = c(-4, 4), link = 'probit', beta = beta[3,], alpha = alpha[2,])
