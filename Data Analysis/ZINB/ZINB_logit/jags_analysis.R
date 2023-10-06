
zipModelkm_logit<- "model{

#beta is the coefficient for log link regression
#alpha is the coefficient for logistic regression

## Likelihood
    for(i in 1:n){
        y[i] ~ dpois(mu[i])
        mu[i] <- lambda[i]*(1-ind[i]) + 1e-10*ind[i]      #mean of zip model

        lambda[i] <- exp(mu.count[i])                     #mean of poisson component
        mu.count[i] <- inprod(beta[],X[i,])             #product of beta vector and X design matrix

        ## Zero-Inflation
        ind[i] ~ dbern(theta[i])                          #indicate where zero count belongs
        theta[i] <- ilogit(mu.binary[i])            #probability zero count is an excess zero
        mu.binary[i] <- inprod(alpha[], Z[i,])      #product of alpha vector and Z design matrix
    }
     
    
    
    
    
    ## Priors for alpha
  
   
    for(j in 1:nKz){
    
    indz[j]~dbern(0.5)   # 0's and 1's
    
    
    pz[j]~dnorm(0,100)
   
   
    
   
        alpha[j] <-(indz[j])*pz[j]
    }
    
    
    
    
     ## Priors for  beta
  
    for(k in 1:nKx){
    
    indx[k]~dbern(0.5)   # 0's and 1's
    
    px[k]~dnorm(0,1/10)
   
   
 
    
   
    
        beta[k] <-(indx[k])*px[k]
    }
    
      
    
    
}"


zipModelkm_clog<-"model{

#beta is the coefficient for log link regression
#alpha is the coefficient for logistic regression

## Likelihood
    for(i in 1:n){
        y[i] ~ dpois(mu[i])
        mu[i] <- lambda[i]*(1-ind[i]) + 1e-10*ind[i]      #mean of zip model

        lambda[i] <- exp(mu.count[i]+log(offset.variable[i]))                     #mean of poisson component
        mu.count[i] <- inprod(beta[],X[i,])             #product of beta vector and X design matrix

        ## Zero-Inflation
        ind[i] ~ dbern(theta[i])                          #indicate where zero count belongs
        log(theta[i]) <- -exp(mu.binary[i])            #probability zero count is an excess zero
        mu.binary[i] <- inprod(alpha[], Z[i,])      #product of alpha vector and Z design matrix
    }
     
    
    
    
    
    ## Priors for alpha
  
   
    for(j in 1:nKz){
    
    indz[j]~dbern(0.5)   # 0's and 1's
    
    
    pz[j]~dnorm(0,100)
   
   
    
   
        alpha[j] <-(indz[j])*pz[j]
    }
    
    
    
    
     ## Priors for  beta
  
    for(k in 1:nKx){
    
    indx[k]~dbern(0.5)   # 0's and 1's
    
    px[k]~dnorm(0,1/1000)
   
        beta[k] <-(indx[k])*px[k]
    }
    
      
    
    
}"




zipModelkm_negbin<-"model{


  
  #beta is the coefficient for log link regression
  #alpha is the coefficient for logistic regression
  
  ## Likelihood
  for(i in 1:n){
        y[i] ~ dnegbin(p[i],size)
        p[i] <- size/(size+(1-ind[i])*lambda[i])- 1e-10*ind[i]  
        lambda[i] <- exp(mu.count[i]+log(offset.variable[i]))                    
        mu.count[i] <- inprod(beta[],X[i,])             #product of beta vector and X design matrix
     
    ## Zero-Inflation
    ind[i] ~ dbern(theta[i])                         #indicate where zero count belongspo
    theta[i] <- ilogit(mu.binary[i])            #probability zero count is an excess zero
    mu.binary[i] <- inprod(alpha[], Z[i,])      #product of alpha vector and Z design matrix
  }
  
  
   
  
    
   ## Priors

   

    for(k.p in 1:Kp){
        beta[k.p] ~ dnorm(0, 0.01)
    }

    for(k.z in 1:Kz){
        alpha[k.z] ~ dnorm(0, 0.01)
    }


    
    for(k in (Kp+1):nKx){
        
        indx[k] ~ dbern(0.5)
        zz[k] <- equals(indx[k], 0)
        F[k] <- zz[k]*0.00001 + 1 - zz[k]
        eta[k] ~ dgamma(0.0001,0.0001)
        gamma[k] <- F[k]*(1/eta[k])
        beta[k] ~ dnorm(0, 1/gamma[k])
    }

    for(j in (Kz+1):nKz){
        
        indz[j] ~ dbern(0.5)
        zza[j] <- equals(indz[j], 0)
        Fa[j] <- zza[j]*0.00001 + 1 - zza[j]
        etaa[j] ~ dgamma(0.0001,0.0001)
        gammaa[j] <- Fa[j]*(1/etaa[j])
        alpha[j] ~ dnorm(0, 1/gammaa[j])
    }

    
  size<-1000
#size ~dunif(0.00001,1000)  



    
}"



zipModelkm_negbintt<-"model{


  
  #beta is the coefficient for log link regression
  #alpha is the coefficient for logistic regression
  
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
        beta[k.p] ~ dnorm(0, 0.01)
    }

    for(k.z in 1:Kz){
        alpha[k.z] ~ dnorm(0, 0.01)
    }


    
    for(k in (Kp+1):nKx){
        
        indx[k] ~ dbern(0.5)
        zz[k] <- equals(indx[k], 0)
        F[k] <- zz[k]*0.00001 + 1 - zz[k]
        eta[k] ~ dgamma(0.0001,0.0001)
        gamma[k] <- F[k]*(1/eta[k])
        beta[k] ~ dnorm(0, 1/gamma[k])
    }

    for(j in (Kz+1):nKz){
        
        indz[j] ~ dbern(0.5)
        zza[j] <- equals(indz[j], 0)
        Fa[j] <- zza[j]*0.00001 + 1 - zza[j]
        etaa[j] ~ dgamma(0.0001,0.0001)
        gammaa[j] <- Fa[j]*(1/etaa[j])
        alpha[j] ~ dnorm(0, 1/gammaa[j])
    }

    
    
#size ~dunif(0.00001,50)  
size<-30


    
}"


