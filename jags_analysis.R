zipModellogit <- "model{
## Likelihood
    for(i in 1:n){
        y[i] ~ dpois(mu[i])
        mu[i] <- lambda[i]*(1-zero[i]) + 1e-10*zero[i]

        lambda[i] <- exp(mu.count[i]+ log(offset.variable[i]))
        mu.count[i] <- inprod(beta[],X[i,]) 
         
        ## Zero-Inflation
        zero[i] ~ dbern(pi[i])
        pi[i] <- ilogit(mu.binary[i])
        mu.binary[i] <- inprod(alpha[], Z[i,])
    }

    ## Priors

   

    for(k.p in 1:Kp){
        beta[k.p] ~ dnorm(0, 0.0001)
    }

    for(k.z in 1:Kz){
        alpha[k.z] ~ dnorm(0, 0.0001)
    }


    
    for(k in (Kp+1):nKx){
        
        indx[k] ~ dbern(0.5)
        zz[k] <- equals(indx[k], 0)
        F[k] <- zz[k]*0.00001 + 1 - zz[k]
        eta[k] ~ dgamma(0.5,0.5)
        gamma[k] <- F[k]*(1/eta[k])
        beta[k] ~ dnorm(0, 1/gamma[k])
    }

    for(j in (Kz+1):nKz){
        
        indz[j] ~ dbern(0.5)
        zza[j] <- equals(indz[j], 0)
        Fa[j] <- zza[j]*0.00001 + 1 - zza[j]
        etaa[j] ~ dgamma(0.5,0.5)
        gammaa[j] <- Fa[j]*(1/etaa[j])
        alpha[j] ~ dnorm(0, 1/gammaa[j])
    }

}"





zipModelkm_clog_off<- "model{

#beta is the coefficient for log link regression
#alpha is the coefficient for logistic regression

## Likelihood
    for(i in 1:n){
        y[i] ~ dpois(mu[i])
        mu[i] <- lambda[i]*(1-ind[i]) + 1e-10*ind[i]      #mean of zip model

       lambda[i] <- exp(mu.count[i]+log(offset.variable[i]))                      #mean of poisson component
        mu.count[i] <- inprod(beta[],X[i,])             #product of beta vector and X design matrix

        ## Zero-Inflation
        ind[i] ~ dbern(theta[i])                          #indicate where zero count belongs
        log(theta[i]) <- -exp(mu.binary[i])            #probability zero count is an excess zero
        mu.binary[i] <- inprod(alpha[], Z[i,])      #product of alpha vector and Z design matrix
    }
     
    
    
    for(k.p in 1:Kp){
        beta[k.p] ~ dnorm(0, 0.0001)
    }

    for(k.z in 1:Kz){
        alpha[k.z] ~ dnorm(0, 0.0001)
    }


    
    for(k in (Kp+1):nKx){
        
        indx[k] ~ dbern(0.5)
        zz[k] <- equals(indx[k], 0)
        F[k] <- zz[k]*0.00001 + 1 - zz[k]
        eta[k] ~ dgamma(0.5,0.5)
        gamma[k] <- F[k]*(1/eta[k])
        beta[k] ~ dnorm(0, 1/gamma[k])
    }

    for(j in (Kz+1):nKz){
        
        indz[j] ~ dbern(0.5)
        zza[j] <- equals(indz[j], 0)
        Fa[j] <- zza[j]*0.00001 + 1 - zza[j]
        etaa[j] ~ dgamma(0.5, 0.5)
        gammaa[j] <- Fa[j]*(1/etaa[j])
        alpha[j] ~ dnorm(0, 1/gammaa[j])
    }
}"


##################################################################
############## ZINB model

ZINBlogit <- "model{
## Likelihood
    ## Likelihood
  for(i in 1:n){
        y[i] ~ dnegbin(p[i],size)
        p[i] <- size/(size+(1-ind[i])*lambda[i])- 1e-10*ind[i]
        lambda[i] <- exp(mu.count[i]+log(offset.variable[i]))                    
        mu.count[i] <- inprod(beta[],X[i,])     
     
    ## Zero-Inflation
    ind[i] ~ dbern(theta[i])                 
    theta[i] <- ilogit(mu.binary[i])            #probability zero count is an excess zero
    mu.binary[i] <- inprod(alpha[], Z[i,])      #product of alpha vector and Z design matrix
  }
  

    ## Priors

   

    for(k.p in 1:Kp){
        beta[k.p] ~ dnorm(0, 0.0001)
    }

    for(k.z in 1:Kz){
        alpha[k.z] ~ dnorm(0, 0.0001)
    }


    
    for(k in (Kp+1):nKx){
        
        indx[k] ~ dbern(0.5)
        zz[k] <- equals(indx[k], 0)
        F[k] <- zz[k]*0.00001 + 1 - zz[k]
        eta[k] ~ dgamma(0.5,0.5)
        gamma[k] <- F[k]*(1/eta[k])
        beta[k] ~ dnorm(0, 1/gamma[k])
    }

    for(j in (Kz+1):nKz){
        
        indz[j] ~ dbern(0.5)
        zza[j] <- equals(indz[j], 0)
        Fa[j] <- zza[j]*0.00001 + 1 - zza[j]
        etaa[j] ~ dgamma(0.5,0.5)
        gammaa[j] <- Fa[j]*(1/etaa[j])
        alpha[j] ~ dnorm(0, 1/gammaa[j])
    }
    
    #size ~dunif(0.00001,50)  
    size<-30

}"



