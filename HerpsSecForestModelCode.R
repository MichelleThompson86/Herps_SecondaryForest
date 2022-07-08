model{
  
  
  ## **********************************************************************
  ## Detection Component
  ## **********************************************************************  
  
  p.mean~dunif(0,1)  # mean detection probability across species
  b  <- logit(p.mean)  # logit of the mean detection probability across species
  
  sigma.v ~ dunif(0,10) #standard deviation (on the logit scale) of the distribution from which the detection probability of each species is drawn
  tau.v <- pow(sigma.v,-2) 
  
  
  ## Effect of survey variables on detection
  mu.alpha1 ~ dnorm(0, 0.01) #across-species mean effect of season on logit detection probability
  mu.alpha2 ~ dnorm(0, 0.01) #across-species mean effect of effort on logit detection probability
  sigma.alpha1 ~ dunif(0,10) # standard deviation (on the logit scale) of the distribution from which the effect of season on detection probability for each species is drawn
  sigma.alpha2 ~ dunif(0,10) # standard deviation (on the logit scale) of the distribution from which the effect of effort on detection probability for each species is drawn
  tau.alpha1 <- 1/pow(sigma.alpha1,2)
  tau.alpha2 <- 1/pow(sigma.alpha2,2)
  
  
  # effect of season and effort on detection (both variable over time)
  for(transect in 1:ntransect) {
    for(socc in 1:nsocc){
      season[transect,socc] ~ dnorm(0, 1) 
      effort[transect,socc] ~ dnorm(0, 1)
    }  #end socc (survey occasion) loop
  }  #end ntransect loop
  
  
  ## species-specific 
  for(sp in 1:nsp) {
    alpha1[sp]~dnorm(mu.alpha1, tau.alpha1) #T(-5,5) # species-specific coefficient for effect of season logit detection probability, can truncate if needed
    alpha2[sp]~dnorm(mu.alpha2, tau.alpha2)#T(-5,5) # species-specific coefficient for effect of standardized effort logit detection probability, can truncate if needed
  } #end sp loop
  
  #random effect of site on detection
  sigma.site.p ~ dunif(0,10)  # if data sparse want more shrinkage can use dexp(1)  or half cauchy
  tau.site.p <- pow(sigma.site.p, -2)
  
  for (site in 1:nsite) {
    eta.site.p[site] ~ dnorm(0,tau.site.p)
  }  # end site loop
  
  
  ## **********************************************************************
  ## Occupancy Component
  ## **********************************************************************  
  
  
  psi.mean ~ dunif(0,1) # mean probability of occurrence across species
  a <- log(psi.mean)-log(1-psi.mean) # logit of the mean probability of occurrence across species
  sigma.u ~ dunif(0,10) # standard deviation (on the logit scale) of the distribution from which the probability of occurrence of each species is drawn
  tau.u <- pow(sigma.u,-2)
  
  
  #random effect of site on occupancy
  sigma.site.psi ~ dunif(0,10) # if data sparse want more shrinkage can use dexp(1)  or half cauchy
  tau.site.psi <- pow(sigma.site.psi, -2)
  
  for(site in 1:nsite){
    eta.site.psi[site] ~ dnorm(0, tau.site.psi) # random intercept for each site
  } # end site loop
  
     
  #Equal vector of probability for tree sampling
  for (k in 1:Ntree) {
    c[k] <- 1/Ntree
  }
  #Tree sampling and variance-covariance matrix construction
  K ~ dcat(c[])         

  
  
  ## incorporate phylogenetic covariance structure 
  
  lambda1 ~ dunif(0,1)
  lambda2 ~ dunif(0,1)  
  lambda3 ~ dunif(0,1)
  lambda4 ~ dunif(0,1)
  lambda5 ~ dunif(0,1)
  lambda6 ~ dunif(0,1)
  lambda7 ~ dunif(0,1)
  lambda8 ~ dunif(0,1)
  lambda9 ~ dunif(0,1)
  
  ## VCOV is an array of variance covariance (or correlation) matrices
  ## ID is a species ID matrix

  
  beta.mat1[1:nsp,1:nsp] <- lambda1*VCOV[,,K] + (1-lambda1)*ID[,]
  beta.mat2[1:nsp,1:nsp] <- lambda2*VCOV[,,K] + (1-lambda2)*ID[,]
  beta.mat3[1:nsp,1:nsp] <- lambda3*VCOV[,,K] + (1-lambda3)*ID[,]
  beta.mat4[1:nsp,1:nsp] <- lambda4*VCOV[,,K] + (1-lambda4)*ID[,]
  beta.mat5[1:nsp,1:nsp] <- lambda5*VCOV[,,K] + (1-lambda5)*ID[,]
  beta.mat6[1:nsp,1:nsp] <- lambda6*VCOV[,,K] + (1-lambda6)*ID[,]
  beta.mat7[1:nsp,1:nsp] <- lambda7*VCOV[,,K] + (1-lambda7)*ID[,]
  beta.mat8[1:nsp,1:nsp] <- lambda8*VCOV[,,K] + (1-lambda8)*ID[,]
  beta.mat9[1:nsp,1:nsp] <- lambda9*VCOV[,,K] + (1-lambda9)*ID[,]
  
  sigma.psi.stage1 ~ dunif(0,10)
  tau.beta.rand.mat1[1:nsp,1:nsp] <-inverse((sigma.psi.stage1^2)*beta.mat1[,])
  sigma.psi.stage2 ~ dunif(0,10)
  tau.beta.rand.mat2[1:nsp,1:nsp] <-inverse((sigma.psi.stage2^2)*beta.mat2[,])
  sigma.psi.stage3 ~ dunif(0,10)
  tau.beta.rand.mat3[1:nsp,1:nsp] <-inverse((sigma.psi.stage3^2)*beta.mat3[,])
  sigma.psi.stage4 ~ dunif(0,10)
  tau.beta.rand.mat4[1:nsp,1:nsp] <-inverse((sigma.psi.stage4^2)*beta.mat4[,])
  sigma.psi.habitat ~ dunif(0,10)
  tau.beta.rand.mat5[1:nsp,1:nsp] <-inverse((sigma.psi.habitat^2)*beta.mat5[,])
  sigma.psi.stage1hab ~ dunif(0,10)
  tau.beta.rand.mat6[1:nsp,1:nsp] <-inverse((sigma.psi.stage1hab^2)*beta.mat6[,])
  sigma.psi.stage2hab ~ dunif(0,10)
  tau.beta.rand.mat7[1:nsp,1:nsp] <-inverse((sigma.psi.stage2hab^2)*beta.mat7[,])
  sigma.psi.stage3hab ~ dunif(0,10)
  tau.beta.rand.mat8[1:nsp,1:nsp] <-inverse((sigma.psi.stage3hab^2)*beta.mat8[,])
  sigma.psi.stage4hab ~ dunif(0,10)
  tau.beta.rand.mat9[1:nsp,1:nsp] <-inverse((sigma.psi.stage4hab^2)*beta.mat9[,])
  
  # Generate host of zeros for MVN  
  for(sp in 1:nsp) {  	
    zeros.sp[sp]<-0
    
  } #end sp loop
  
  ## draw species-specific slopes
  psi.beta.resid1[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat1[,])
  ## draw species-specific slopes
  psi.beta.resid2[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat2[,])
  ## draw species-specific slopes
  psi.beta.resid3[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat3[,])
  ## draw species-specific slopes
  psi.beta.resid4[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat4[,])
  ## draw species-specific slopes
  psi.beta.resid5[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat5[,])
  ## draw species-specific slopes
  psi.beta.resid6[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat6[,])
  ## draw species-specific slopes
  psi.beta.resid7[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat7[,])
  ## draw species-specific slopes
  psi.beta.resid8[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat8[,])
  ## draw species-specific slopes
  psi.beta.resid9[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat9[,])
  
  # Mean response to environment across all species
  mu.beta1  ~ dnorm(0, 0.01)
  mu.beta2  ~ dnorm(0, 0.01)
  mu.beta3  ~ dnorm(0, 0.01)
  mu.beta4  ~ dnorm(0, 0.01)
  mu.beta5  ~ dnorm(0, 0.01)
  mu.beta6  ~ dnorm(0, 0.01)
  mu.beta7  ~ dnorm(0, 0.01)
  mu.beta8  ~ dnorm(0, 0.01)
  mu.beta9  ~ dnorm(0, 0.01)
  
  # Create total psi.beta.sp by adding phylogenetic residual to mean response
  for (sp in 1:nsp){
    psi.beta.sp1[sp] <- mu.beta1 + psi.beta.resid1[sp] 
    psi.beta.sp2[sp] <- mu.beta2 + psi.beta.resid2[sp]
    psi.beta.sp3[sp] <- mu.beta3 + psi.beta.resid3[sp]
    psi.beta.sp4[sp] <- mu.beta4 + psi.beta.resid4[sp]
    psi.beta.sp5[sp] <- mu.beta5 + psi.beta.resid5[sp]
    psi.beta.sp6[sp] <- mu.beta6 + psi.beta.resid6[sp]
    psi.beta.sp7[sp] <- mu.beta7 + psi.beta.resid7[sp]
    psi.beta.sp8[sp] <- mu.beta8 + psi.beta.resid8[sp]
    psi.beta.sp9[sp] <- mu.beta9 + psi.beta.resid9[sp]
  }
  
  #cov det-occ  
  rho~dunif(-1,1) # a measure of the degree of covariance between probability of occurrence and detection probability (thought to be mediated by abundance)
  var.v<-tau.v/(1-pow(rho,2)) # adjustment of precision of detection probability to account for covariance between detection and occurrence
  
  for(sp in 1:nsp) {
    u[sp] ~ dnorm(a,tau.u)T(-5,5) # species-specific mean probability of occurrence
    mu.v[sp] <- b + (rho*sigma.v/sigma.u)*(u[sp]-a) # sets up the correlation between probability of occurrence and detection probability
    v[sp] ~ dnorm(mu.v[sp],var.v)T(-10,5) # species-specific mean detection probability
    
    
    #estimate the occupancy probability 
    for(transect in 1:ntransect) {
      logit(psi[transect,sp]) <-
        u[sp]+                           #species intercept
        psi.beta.sp1[sp]*stage1[transect] + psi.beta.sp2[sp]*stage2[transect] + psi.beta.sp3[sp]*stage3[transect] + psi.beta.sp4[sp]*stage4[transect] + 
        psi.beta.sp5[sp]*habitat[transect]+ psi.beta.sp6[sp]*habitat[transect]*stage1[transect]+ psi.beta.sp7[sp]*habitat[transect]*stage2[transect]+
        psi.beta.sp8[sp]*habitat[transect]*stage3[transect]+psi.beta.sp9[sp]*habitat[transect]*stage4[transect] +
        eta.site.psi[site[transect]]
      
      # Generate occupancy states
      mu.psi[transect,sp]<-psi[transect,sp]
      Z[transect,sp] ~ dbern(mu.psi[transect,sp])
      
      for(socc in 1:nsocc){
        logit(p[transect,socc,sp]) <- 
          v[sp] + 
          alpha1[sp]*season[transect, socc] +
          alpha2[sp]* effort[transect,socc] + 
          eta.site.p[site[transect]]
        
        
        # Generate expectation of detection
        mu.p[transect, socc, sp] <- Z[transect, sp]* p[transect, socc, sp] 
        
        # Compare expectation to data
        X[transect,socc,sp] ~ dbern(mu.p[transect,socc,sp])
        
      } #end socc loop
    } #end transect loop
  } #end sp loop
  
  # derived parameters
  for(transect in 1:ntransect){
    SpRT[transect] <- sum(Z[transect,]) # determine species richness at each transect by summing across species indicators
  }
  
} # model.jags


#############################with trait##############################

model{
  
  
  ## **********************************************************************
  ## Detection Component
  ## **********************************************************************  
  
  p.mean~dunif(0,1)  # mean detection probability across species
  b  <- logit(p.mean)  # logit of the mean detection probability across species
  
  sigma.v ~ dunif(0,10) #standard deviation (on the logit scale) of the distribution from which the detection probability of each species is drawn
  tau.v <- pow(sigma.v,-2) 
  
  
  ## Effect of survey variables on detection
  mu.alpha1 ~ dnorm(0, 0.01) #across-species mean effect of season on logit detection probability
  mu.alpha2 ~ dnorm(0, 0.01) #across-species mean effect of effort on logit detection probability
  sigma.alpha1 ~ dunif(0,10) # standard deviation (on the logit scale) of the distribution from which the effect of season on detection probability for each species is drawn
  sigma.alpha2 ~ dunif(0,10) # standard deviation (on the logit scale) of the distribution from which the effect of effort on detection probability for each species is drawn
  tau.alpha1 <- 1/pow(sigma.alpha1,2)
  tau.alpha2 <- 1/pow(sigma.alpha2,2)
  
  
  # effect of season and effort on detection (both variable over time)
  for(transect in 1:ntransect) {
    for(socc in 1:nsocc){
      season[transect,socc] ~ dnorm(0, 1) 
      effort[transect,socc] ~ dnorm(0, 1)
    }
  }
  
  
  ## species-specific 
  for(sp in 1:nsp) {
    alpha1[sp]~dnorm(mu.alpha1, tau.alpha1) #T(-5,5) # species-specific coefficient for effect of season logit detection probability
    alpha2[sp]~dnorm(mu.alpha2, tau.alpha2) #T(-5,5) # species-specific coefficient for effect of standardized effort logit detection probability
  }
  
  #random effect of site on detection
  sigma.site.p ~ dunif(0,10) 
  tau.site.p <- pow(sigma.site.p, -2)
  
  for (site in 1:nsite) {
    eta.site.p[site] ~ dnorm(0,tau.site.p)
  }
  
  
  ## **********************************************************************
  ## Occupancy Component
  ## **********************************************************************  
  
  
  psi.mean ~ dunif(0,1) # mean probability of occurrence across species
  a <- log(psi.mean)-log(1-psi.mean) # logit of the mean probability of occurrence across species
  sigma.u ~ dunif(0,10) # standard deviation (on the logit scale) of the distribution from which the probability of occurrence of each species is drawn
  tau.u <- pow(sigma.u,-2)
  
  
  #random effect of site on occupancy
  sigma.site.psi ~ dunif(0,10)
  tau.site.psi <- pow(sigma.site.psi, -2)
  for(site in 1:nsite){
    eta.site.psi[site] ~ dnorm(0, tau.site.psi) # random intercept for each site
  }
  
        
  #Equal vector of probability for tree sampling
  for (k in 1:Ntree) {
    c[k] <- 1/Ntree
  }
  #Tree sampling and variance-covariance matrix construction
  K ~ dcat(c[])         

  
  ## incorporate phylogenetic covariance structure
  
  lambda1 ~ dunif(0,1)
  lambda2 ~ dunif(0,1)  
  lambda3 ~ dunif(0,1)
  lambda4 ~ dunif(0,1)
  lambda5 ~ dunif(0,1)
  lambda6 ~ dunif(0,1)
  lambda7 ~ dunif(0,1)
  lambda8 ~ dunif(0,1)
  lambda9 ~ dunif(0,1)
  
  ## VCOV is an array of variance covariance (or correlation) matrices, 
   ## ID is a species ID matrix
  
  beta.mat1[1:nsp,1:nsp] <- lambda1*VCOV[,,K] + (1-lambda1)*ID[,]
  beta.mat2[1:nsp,1:nsp] <- lambda2*VCOV[,,K] + (1-lambda2)*ID[,]
  beta.mat3[1:nsp,1:nsp] <- lambda3*VCOV[,,K] + (1-lambda3)*ID[,]
  beta.mat4[1:nsp,1:nsp] <- lambda4*VCOV[,,K] + (1-lambda4)*ID[,]
  beta.mat5[1:nsp,1:nsp] <- lambda5*VCOV[,,K] + (1-lambda5)*ID[,]
  beta.mat6[1:nsp,1:nsp] <- lambda6*VCOV[,,K] + (1-lambda6)*ID[,]
  beta.mat7[1:nsp,1:nsp] <- lambda7*VCOV[,,K] + (1-lambda7)*ID[,]
  beta.mat8[1:nsp,1:nsp] <- lambda8*VCOV[,,K] + (1-lambda8)*ID[,]
  beta.mat9[1:nsp,1:nsp] <- lambda9*VCOV[,,K] + (1-lambda9)*ID[,]
  
  sigma.psi.stage1 ~ dunif(0,10)
  tau.beta.rand.mat1[1:nsp,1:nsp] <-inverse((sigma.psi.stage1^2)*beta.mat1[,])
  sigma.psi.stage2 ~ dunif(0,10)
  tau.beta.rand.mat2[1:nsp,1:nsp] <-inverse((sigma.psi.stage2^2)*beta.mat2[,])
  sigma.psi.stage3 ~ dunif(0,10)
  tau.beta.rand.mat3[1:nsp,1:nsp] <-inverse((sigma.psi.stage3^2)*beta.mat3[,])
  sigma.psi.stage4 ~ dunif(0,10)
  tau.beta.rand.mat4[1:nsp,1:nsp] <-inverse((sigma.psi.stage4^2)*beta.mat4[,])
  sigma.psi.habitat ~ dunif(0,10)
  tau.beta.rand.mat5[1:nsp,1:nsp] <-inverse((sigma.psi.habitat^2)*beta.mat5[,])
  sigma.psi.stage1hab ~ dunif(0,10)
  tau.beta.rand.mat6[1:nsp,1:nsp] <-inverse((sigma.psi.stage1hab^2)*beta.mat6[,])
  sigma.psi.stage2hab ~ dunif(0,10)
  tau.beta.rand.mat7[1:nsp,1:nsp] <-inverse((sigma.psi.stage2hab^2)*beta.mat7[,])
  sigma.psi.stage3hab ~ dunif(0,10)
  tau.beta.rand.mat8[1:nsp,1:nsp] <-inverse((sigma.psi.stage3hab^2)*beta.mat8[,])
  sigma.psi.stage4hab ~ dunif(0,10)
  tau.beta.rand.mat9[1:nsp,1:nsp] <-inverse((sigma.psi.stage4hab^2)*beta.mat9[,])
  
  # Generate host of zeros for MVN  
  for(sp in 1:nsp) {  	
    zeros.sp[sp]<-0
  } 
  
  ## draw species-specific slopes
  psi.beta.resid1[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat1[,])
  ## draw species-specific slopes
  psi.beta.resid2[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat2[,])
  ## draw species-specific slopes
  psi.beta.resid3[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat3[,])
  ## draw species-specific slopes
  psi.beta.resid4[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat4[,])
  ## draw species-specific slopes
  psi.beta.resid5[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat5[,])
  ## draw species-specific slopes
  psi.beta.resid6[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat6[,])
  ## draw species-specific slopes
  psi.beta.resid7[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat7[,])
  ## draw species-specific slopes
  psi.beta.resid8[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat8[,])
  ## draw species-specific slopes
  psi.beta.resid9[1:nsp] ~ dmnorm(zeros.sp[], tau.beta.rand.mat9[,])
  
  # Mean response to environment across all species
  mu.beta1  ~ dnorm(0, 0.01)
  mu.beta2  ~ dnorm(0, 0.01)
  mu.beta3  ~ dnorm(0, 0.01)
  mu.beta4  ~ dnorm(0, 0.01)
  mu.beta5  ~ dnorm(0, 0.01)
  mu.beta6  ~ dnorm(0, 0.01)
  mu.beta7  ~ dnorm(0, 0.01)
  mu.beta8  ~ dnorm(0, 0.01)
  mu.beta9  ~ dnorm(0, 0.01)
  psi.beta.trait1 ~ dnorm(0, 0.01)
  psi.beta.trait2 ~ dnorm(0, 0.01)
  psi.beta.trait3 ~ dnorm(0, 0.01)
  psi.beta.trait4 ~ dnorm(0, 0.01)
  psi.beta.trait5 ~ dnorm(0, 0.01)
  psi.beta.trait6 ~ dnorm(0, 0.01)
  psi.beta.trait7 ~ dnorm(0, 0.01)
  psi.beta.trait8 ~ dnorm(0, 0.01)
  psi.beta.trait9 ~ dnorm(0, 0.01)
  
  # Create total psi.beta.sp by adding phylogenetic residual to mean response
  for (sp in 1:nsp){
    psi.beta.sp1[sp] <- mu.beta1 + psi.beta.trait1*trait[sp] + psi.beta.resid1[sp] 
    psi.beta.sp2[sp] <- mu.beta2 + psi.beta.trait2*trait[sp] + psi.beta.resid2[sp] 
    psi.beta.sp3[sp] <- mu.beta3 + psi.beta.trait3*trait[sp] + psi.beta.resid3[sp] 
    psi.beta.sp4[sp] <- mu.beta4 + psi.beta.trait4*trait[sp] + psi.beta.resid4[sp] 
    psi.beta.sp5[sp] <- mu.beta5 + psi.beta.trait5*trait[sp] + psi.beta.resid5[sp] 
    psi.beta.sp6[sp] <- mu.beta6 + psi.beta.trait6*trait[sp] + psi.beta.resid6[sp] 
    psi.beta.sp7[sp] <- mu.beta7 + psi.beta.trait7*trait[sp] + psi.beta.resid7[sp] 
    psi.beta.sp8[sp] <- mu.beta8 + psi.beta.trait8*trait[sp] + psi.beta.resid8[sp] 
    psi.beta.sp9[sp] <- mu.beta9 + psi.beta.trait9*trait[sp] + psi.beta.resid9[sp] 
  }
  
  #cov det-occ  
  rho~dunif(-1,1) # a measure of the degree of covariance between probability of occurrence and detection probability (thought to be mediated by abundance)
  var.v<-tau.v/(1-pow(rho,2)) # adjustment of precision of detection probability to account for covariance between detection and occurrence
  
  for(sp in 1:nsp) {
    u[sp] ~ dnorm(a,tau.u)T(-5,5) # species-specific mean probability of occurrence
    mu.v[sp] <- b + (rho*sigma.v/sigma.u)*(u[sp]-a) # sets up the correlation between probability of occurrence and detection probability
    v[sp] ~ dnorm(mu.v[sp],var.v)T(-10,5) # species-specific mean detection probability
    
    
    #estimate the occupancy probability 
    for(transect in 1:ntransect) {
      logit(psi[transect,sp]) <-
        u[sp]+                           #species intercept
        psi.beta.sp1[sp]*stage1[transect] + psi.beta.sp2[sp]*stage2[transect] + psi.beta.sp3[sp]*stage3[transect] + psi.beta.sp4[sp]*stage4[transect] + 
        psi.beta.sp5[sp]*habitat[transect]+ psi.beta.sp6[sp]*habitat[transect]*stage1[transect]+ psi.beta.sp7[sp]*habitat[transect]*stage2[transect]+
        psi.beta.sp8[sp]*habitat[transect]*stage3[transect]+psi.beta.sp9[sp]*habitat[transect]*stage4[transect] +
        eta.site.psi[site[transect]]
      
      # Generate occupancy states
      mu.psi[transect,sp]<-psi[transect,sp]
      Z[transect,sp] ~ dbern(mu.psi[transect,sp])
      
      for(socc in 1:nsocc){
        logit(p[transect,socc,sp]) <- 
          v[sp] + 
          alpha1[sp]*season[transect, socc] +
          alpha2[sp]* effort[transect,socc] + 
          eta.site.p[site[transect]]
        
        
        
        # Generate expectation of detection
        mu.p[transect, socc, sp] <- Z[transect, sp]* p[transect, socc, sp] 
        
        # Compare expectation to data
        X[transect,socc,sp] ~ dbern(mu.p[transect,socc,sp])
        
      }#/socc
    } #/transect
  } #/sp
  
  
} # /model.jags

