zipMod<- "model{

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
     
    
 #priors for alpha
 
 
  
  
 for(k in 1:nKz){
 
 indz[k]~dbern(0.5)   # 0's and 1's
    
    fz[k]<-indz[k]*0.01+1-indz[k]
    
    #if indz[k]=1 then fz[k]=0.001 which is a spike
    #if indz[k]=0 then fz[k]=1 which is a slab
   
   
    tauz[k]~dgamma(10E-3,10E-3)
    
    #sigma=1/tau
    
    varz[k]<-fz[k]*(1/tauz[k])
    
        alpha[k] ~ dnorm(0,varz[k])
 }
 
 
 
 
 #priors for beta
 
 
 
 for(j in 1:nKx){
 indx[j]~dbern(0.5)   # 0's and 1's
    
    fx[j]<-indx[j]*0.01+1-indx[j]
    
    #if indx[j]=1 then fx[j]=0.01 which is a spike
    #if indx[j]=0 then fx[j]=1 which is a slab
   
   
    taux[j]~dgamma(10E-3,10E-3)
    
    #sigma=1/tau
    
    varx[j]<-fx[j]*(1/taux[j])
    
       beta[j] ~ dnorm(0,varx[j])
 
 }
 
    
}"


zipModel_spikeslab<- "model{

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
    
    indz[j]~dbern(0.4)   
    
    
    pz[j]~dnorm(0,1/100)
    
        alpha[j] <-(indz[j])*pz[j]
        
    }
    
    
    
    
     ## Priors for  beta
  
    for(k in 1:nKx){
    
    indx[k]~dbern(0.4)   
    
    px[k]~dnorm(0,1/100)
   
        beta[k] <-(indx[k])*px[k]
    }
    
      
    
    
}"


zipModel_horseshoe<- "model{

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
  
    #tauz ~ dt(0, 1, 1)T(0,)
    #taux ~ dt(0, 1, 1)T(0,)
    for(j in 1:nKz){
     lz[j] ~ dt(0, 1, 1)T(0,)
     tauz[j] ~ dt(0, 1/nKz, 1)T(0,)
     sigmaz[j] <-  (lz[j])*tauz[j]
     alpha[j] ~ dnorm(0, 1/(sigmaz[j]))
     
     
    }
    
    
    
    
     ## Priors for  beta
  
    for(k in 1:nKx){
     lx[k] ~ dt(0, 1, 1)T(0,)
     taux[k] ~ dt(0, 1/nKx, 1)T(0,)
     sigmax[k] <- (lx[k])*taux[k]
     beta[k] ~ dnorm(0, 1/(sigmax[k]))
    }

}"



