# Bradke. D. R., B. A. Crawford, M. Kaylor, and J. C. Maerz. 2023. Facing uncertainty: evaluating and improving a standard monitoring method for diamond-backed terrapins. Journal of Wildlife Management.

# The "SIMULATION CODE" is for simulating diamond-backed terrapin capture-recapture datasets with a positive change in survival. 
# The "MODEL CODE" is used to evaluate if the simulated increases in survival would be detected using standard monitoring methods.

# load required packages
library(magrittr)
library(nimble)
library(coda)
library(parallel)
library(MCMCvis)

#Read in data files
cov <- read.csv("Tide_cov_2010_2022.csv") #Actual tide data 
MCMCsamp <- read.csv("Mod6_fixed_only_samples_for_sim.csv") #MCMC samples from Jekyll RD CJS model

#############################################
########### SIMULATION CODE##################
#############################################

set.seed(9661)

iter <-100  # number of iterations of the simulation/model to run

### Define simulation settings ###

N1 <-15  # Number of primary capture occasions
per.prime <- 3  # Number of secondary capture occasions in each primary period
N2 <- N1 * per.prime      # Total number of secondary capture occasions
nss <- rep(per.prime,N1) # nss[i] = number of secondary periods in primary period i
cnss <- c(0,cumsum(nss)) %>% head(-1) # cnss[i] = cumulative sum of nss from 1 to i-1

pre.phi <- 0.57  #pre-regulation phi
post.phi <- 0.77 #post-regulation phi

pre.int <- 7  #number of pre-regulation sampling intervals
post.int <- 7  #number of post-regulation sampling intervals

NS <- 60 #Number of individuals existing at start
NN <- round((1-pre.phi)*NS) #Number of new individuals per year (replacement rate = 1-mean survival pre-reg for all scenarios)
NT <- NS + (NN*(N1-1)) #Total number of individuals "existing" for all years

#to create arrays for storing simulation parameters and model results
row <- array(dim=c(iter))
p.mean <- array(dim=c(iter))
p.sd <- array(dim=c(iter))
c.beta <- array(dim=c(iter))
gam.mean <- array(dim=c(iter))
gam.beta <- array(dim=c(iter))
phi.sd <- array(dim=c(iter))

phi.beta.mean <- array(dim=c(iter))
phi.beta.q2.5 <- array(dim=c(iter))
phi.beta.q97.5 <- array(dim=c(iter))

phi.pre.mean <- array(dim=c(iter))
phi.pre.q2.5 <- array(dim=c(iter))
phi.pre.q97.5 <- array(dim=c(iter))

phi.post.mean <- array(dim=c(iter))
phi.post.q2.5 <- array(dim=c(iter))
phi.post.q97.5 <- array(dim=c(iter))

phi.diff.mean <- array(dim=c(iter))
phi.diff.q2.5 <- array(dim=c(iter))
phi.diff.q97.5 <- array(dim=c(iter))

#to draw random MCMC iterations and associated parameter values from the RD CJS model results to be used in simulation iterations
for (k in 1:iter){
  row[k] <- sample(row(MCMCsamp),1,replace=F)
  p.mean[k] <- MCMCsamp$p.mean[row[k]]
  p.sd[k] <- MCMCsamp$p.sd[row[k]]
  c.beta[k] <- MCMCsamp$delta[row[k]]
  gam.mean[k] <- MCMCsamp$gam.mean[row[k]]
  gam.beta[k] <- MCMCsamp$beta[row[k]]
  phi.sd[k] <- MCMCsamp$phi.sd[row[k]]
} #k

#To save simulation settings
Simulation <- cbind(row,p.mean,p.sd,c.beta,gam.mean,gam.beta,N1,N2,NS,NN,pre.phi,post.phi,phi.sd, pre.int,post.int, deparse.level = 1)

write.csv(Simulation,"Simulation_params.csv")

#create file names to store simulation results later
file.num <- c(1:iter)

file.samps <- paste("./Samples/Simulation_samples",file.num,".csv",sep = "")
file.summs <- paste("./Summary/Simulation_summary",file.num,".csv",sep = "")

start.time <- Sys.time()

for (k in 1:iter){
  logit.phi <- c(rnorm(pre.int,qlogis(pre.phi),phi.sd[k]),rnorm(post.int,qlogis(post.phi),phi.sd[k]))
  phi <- plogis(logit.phi) #phi for each sampling interval
  
  logit.p <- rnorm(N2,qlogis(p.mean[k]),p.sd[k])
  
  tide.actual <- cov$Actual_extrap_Jekyll #bring in actual tide vals from sampling at Jekyll during years 1-10 and 12-13 (no sampling in year 11)
  tide.mean <- rep(NA,12)
  cnss.actual <- c(0,3,6,9,12,15,18,21,24,27,30,32,34)
  for(i in 1:10){
    tide.mean[i] <- mean(tide.actual[(cnss.actual[i]+1):(cnss.actual[i]+3)])
  }
  for(i in 12:13){
    tide.mean[i-1] <- mean(tide.actual[(cnss.actual[i]+1):(cnss.actual[i]+2)])
  }
  
  tide <- rnorm((N1-1), mean(tide.mean),sd(tide.mean)) #select tide covariate vals for simulation
  
  logit.gam <- qlogis(gam.mean[k])+gam.beta[k]*tide
  gam <- plogis(logit.gam)
  
  # Define matrices with survival and capture probabilities
  PHI <- matrix(phi, ncol = N1-1, nrow = NT, byrow = TRUE)
  LOGIT.P <- matrix(logit.p, ncol = N2, nrow = NT, byrow = TRUE)
  GAM <- matrix(gam, ncol = N1-1, nrow = NT, byrow = TRUE)
  C.BETA <- c.beta[k]
  
  # Define function to simulate a capture-history (CH) matrix
  simul.rd <- function(PHI, LOGIT.P, C.BETA, GAM, NT, nss, cnss){
    N1 <- dim(PHI)[2] + 1
    N2 <- dim(LOGIT.P)[2]
    
    CH <- matrix(0, ncol = N2, nrow = NT)
    exist.occ <- c(rep(1,NS), rep(2:N1,each = NN)) #vector of occasions when individuals first "exist"
    
    # Fill the CH matrix
    for (i in 1:NT){
      
      #simulate if captured in first secondary occasion 
      #of first primary period that individual exists
      rp <- rbinom(1, 1, plogis(LOGIT.P[i,(cnss[exist.occ[i]]+1)]))
      if (rp==1) {CH[i,(cnss[exist.occ[i]]+1)] <- 1}
      #simulate if captured in remaining secondary occasions 
      #of first primary period that individual exists
      for (j in 2:nss[exist.occ[i]]){
        P <- plogis(LOGIT.P[i,(cnss[exist.occ[i]]+j)]+C.BETA*max(CH[i,(cnss[exist.occ[i]]+1):(cnss[exist.occ[i]]+j-1)]))
        rp <- rbinom(1, 1, P)
        if (rp==1) {CH[i,(cnss[exist.occ[i]]+j)] <- 1}
      }#j
      
      if (exist.occ[i]==N1) {next} #If first exists in final primary period, skip to next individual
      for (t in (exist.occ[i]+1):N1){
        sur <- rbinom(1, 1, PHI[i,(t-1)]) #did individual survive to second primary period?
        if (sur==0) {break}	# If dead, move to next individual 	
        fid <- rbinom(1, 1, GAM[i,(t-1)]) #is individual available in second primary period?
        if (fid==0) {next}
        #if available, was individual captured in first secondary occasion?
        rp <- rbinom(1, 1, plogis(LOGIT.P[i,(cnss[t]+1)]))
        if (rp==1) {CH[i,(cnss[t]+1)] <- 1}         
        for (j in 2:nss[t]){
          #Was individual captured in each remaining secondary occasion? 
          P <- plogis(LOGIT.P[i,(cnss[t]+j)]+C.BETA*max(CH[i,(cnss[t]+1):(cnss[t]+j-1)]))
          rp <- rbinom(1, 1, P)
          if (rp==1) {CH[i,(cnss[t]+j)] <- 1}
        }#j
      } #t
    } #i
    return(CH)
  }
  
  # Execute function
  CH <- simul.rd(PHI, LOGIT.P, C.BETA, GAM, NT, nss, cnss) 
  
  CH <- CH[rowSums(CH)>0, ] #delete individuals with all "0" capture histories
  
  
  #############################################
  ################ MODEL CODE##################
  #############################################
  
  #To get matrix of capture-recapture data for primary periods only
  z.data <- matrix(0, ncol = N1, nrow = dim(CH)[1])
  for (i in 1:dim(CH)[1]){
    for (t in 1:N1){
      cap <- sum(CH[i,(cnss[t]+1):(cnss[t]+per.prime)])
      if (cap>0) z.data[i,t] <- 1
    } #t
  } #i
  
  #To get "if captured before" per primary period (covariate for re-capture effect)
  X <- array(0,dim=dim(CH))
  for (i in 1:dim(CH)[1]){
    for (t in 1:N1){
      for (j in 2:nss[t]) { 
        if(CH[i,(cnss[t]+j-1)]==1) X[i,(cnss[t]+j)] <- 1
        if(X[i,(cnss[t]+j-1)]==1) X[i,(cnss[t]+j)] <- 1
      }
    }
  }
  
  #To get occasion of first capture for each individual
  get.first <- function(x) min(which(x!=0))
  f <- apply(z.data, 1, get.first)
  
  #function for known inclusion states based on capture histories
  known.state.cjs <- function(caphist){
    state <- caphist
    for (i in 1:dim(caphist)[1]){
      n1 <- min(which(caphist[i,]==1))
      n2 <- max(which(caphist[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
    }
    state[state==0] <- NA
    return(state)
  }
  
  #known inclusion states based on capture histories
  z.known <- known.state.cjs(z.data)
  
  #function of initial values for z states
  cjs.init.z <- function(caphist,f){
    for (i in 1:dim(caphist)[1]){
      if (sum(caphist[i,])==1) next
      n2 <- max(which(caphist[i,]==1))
      caphist[i,f[i]:n2] <- NA
    }
    for (i in 1:dim(caphist)[1]){ 
      caphist[i,1:f[i]] <- NA
    }
    return(caphist)
  }
  # initial values for z state
  z.inits <- cjs.init.z(z.data,f)
  
  # known states for availability based on capture histories
  s.known <- as.matrix(z.data)
  s.known[s.known==0]<-NA
  
  BRD <- c(rep(0,pre.int),rep(1,post.int)) #Create binary temporal co-variate indicating pre (0) vs post(1) regulations
  
  
  mod6_code <- nimbleCode({
    # RD model-------------------------------------
    #15 primary occasions, each with 3 secondary sampling periods
    
    #Model 6 - Gamma tide, p rand day, p/c
    # Survival: random year effect and BRD effect
    # Availability (gamma): fixed effect of tide amplitude
    # Capture: random day effect
    # Recapture: not equal to capture
    
    #-------------
    #Section 1. Define priors for all parameters#
    #-------------
    
    #Survival parameters
    for (i in 1:n) {
      for (t in f[i]:(n1-1)){
        logit(phi[i,t]) <- phi.mu + phi.beta * BRD[t] + phi.eps[t] 
      } #t
    } #i
    
    for (t in 1:(n1-1)){
      phi.eps[t] ~ dnorm(0,phi.tau)
    } #t
    
    phi.beta ~ dnorm(0, 0.37)     #Prior for BRD effect
    phi.mean  ~ dunif(0, 1)       #Prior for mean phi
    phi.mu <- log(phi.mean/(1-phi.mean)) 
    phi.tau <- pow(phi.sd,-2)
    phi.sd ~ dunif(0,5) #Prior for sd (random year effect) on survival 
    
    #Availability
    
    for (i in 1:n) {
      for (t in f[i]:(n1-1)){
        logit(gamma[i,t]) <- gam.mu + gam.beta*tide.mean[t]
      } #t
    } #i
    
    gam.mean ~ dunif(0,1) #Prior for mean availability
    gam.mu <- log(gam.mean/(1-gam.mean)) 
    gam.beta ~ dnorm(0,0.37) #Prior for fixed tide effect on availability
    
    #Capture parameters
    
    for (j in 1:n2){
      p.eps[j] ~ dnorm(0,p.tau)
    } #j
    
    for (i in 1:n) {
      for (j in 1:n2) {
        logit(p[i,j]) <- p.mu + p.eps[j] + c.beta * X[i,j] 
      } #j
    } #i
    
    p.mean ~ dunif(0,1) #Prior for mean detection
    p.mu <- log(p.mean/(1-p.mean)) 
    p.tau <- pow(p.sd,-2)
    p.sd ~ dunif(0,5) #Prior for sd (random day effect) on capture probability 
    c.beta ~ dnorm(0, 0.37) #Prior for within secondary period recapture effect (i.e., trap shy effect)
    
    #-------------
    #Section 2. Likelihoods#
    #-------------
    
    # State process #
    for (i in 1:n){
      
      z[i,f[i]] <- 1
      
      
      # Observation process for secondary periods in first primary period of capture
      for (j in 1:nss[f[i]]) { #Loop over secondary periods (times / days)
        y[i,(cnss[f[i]]+j)] ~ dbern(p.eff[i,(cnss[f[i]]+j)])
        p.eff[i,(cnss[f[i]]+j)] <- p[i,(cnss[f[i]]+j)]
      } #j secondary
      
      # State and observation process in subsequent primary periods
      for (t in (f[i]+1):n1){
        phi.eff[i,t] <- phi[i, t-1] * z[i,t-1]  # Effective survival rate
        z[i,t] ~ dbern(phi.eff[i,t])
        avail[i,t] <- gamma[i,t-1] * z[i,t] 
        s[i,t] ~ dbern(avail[i,t])
        
        for (j in 1:nss[t]) { #Loop over secondary periods (times / days)
          y[i,(cnss[t]+j)] ~ dbern(p.eff[i,(cnss[t]+j)])
          p.eff[i,(cnss[t]+j)] <- s[i,t] * p[i,(cnss[t]+j)]
          
        } #j secondary
      } #t primary
    } #i individual
    
    #---------
    #Section 3. Derived parameters#
    #---------
    
    #phi pre-reg, post-reg, and difference between the two
    phi.pre <- 1/(1+exp(-(phi.mu + phi.beta * 0)))
    phi.post <- 1/(1+exp(-(phi.mu + phi.beta * 1)))
    phi.diff <- phi.post-phi.pre #effect size of BRD reg on survival
    
  }) #End Model
  
  #constants
  Mod6Consts <- list(f = f, n=nrow(CH),n1=N1,n2=N2,nss=nss,cnss=cnss)
  #data
  Mod6Data <- list(y=CH,tide.mean=tide,z=z.known, BRD=BRD, s=s.known,X=X)
  
  #Parameters monitored
  params<- c("phi.mean","phi.sd","phi.beta","phi.pre","phi.post","phi.diff","p.mean","p.sd","c.beta","gam.beta","gam.mean")
  
  #Set initial values
  inits <- list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),phi.beta=rnorm(1),gam.mean=runif(1,0,1),gam.beta=rnorm(1),p.mean=runif(1,0,1),c.beta=rnorm(1),p.sd=runif(1,0,5),z=z.inits)
  
  #Run in parallel
  cl <- makeCluster(3)
  clusterExport(cl = cl, varlist = c("Mod6Consts", "Mod6Data", "inits", "params", "mod6_code"))
  mod6.out <- clusterEvalQ(cl = cl,{
    library(nimble) #you're now in a totally different environment so have to load the package again
    library(coda)
    Model6 <- nimbleModel(code = mod6_code, name = 'mod6', constants = Mod6Consts, data = Mod6Data, inits = inits)
    Model6$initializeInfo()
    mcmcMod6 <- configureMCMC(Model6, monitors = params, print = T)
    Mod6MCMC <- buildMCMC(mcmcMod6) #actually build the code for those samplers
    Cmodel <- compileNimble(Model6) #compiling the model itself in C++;
    CMod6MCMC <- compileNimble(Mod6MCMC, project = Model6) # compile the samplers next
    CMod6MCMC$run(niter = 120000, nburnin = 30000)
    
    return(as.mcmc(as.matrix(CMod6MCMC$mvSamples)))
    
  })
  
  
  Samp.mod6 <- mcmc.list(mod6.out)
  
  stopCluster(cl)
  
  #Save summary and and posterior parameter estimates (MCMC samples)
  
  summary <- MCMCsummary(mod6.out)
  write.csv(summary, file=file.summs[k])
  
  samples <- as.data.frame(rbind(mod6.out[[1]],mod6.out[[2]],mod6.out[[3]]))
  write.csv(samples, file=file.samps[k])
  

  #compile results  
  
  phi.beta.mean[k] <- summary['phi.beta',1]
  phi.beta.q2.5[k] <- summary['phi.beta',3]
  phi.beta.q97.5[k] <- summary['phi.beta',5]
  
  phi.pre.mean[k] <- summary['phi.pre',1]
  phi.pre.q2.5[k] <- summary['phi.pre',3]
  phi.pre.q97.5[k] <- summary['phi.pre',5]
  
  phi.post.mean[k] <- summary['phi.post',1]
  phi.post.q2.5[k] <- summary['phi.post',3]
  phi.post.q97.5[k] <- summary['phi.post',5]
  
  phi.diff.mean[k] <- summary['phi.diff',1]
  phi.diff.q2.5[k] <- summary['phi.diff',3]
  phi.diff.q97.5[k] <- summary['phi.diff',5]
  
}

Sys.time() - start.time

Simulation.results <- cbind(phi.pre.mean, phi.pre.q2.5, phi.pre.q97.5, phi.post.mean, phi.post.q2.5,phi.post.q97.5, phi.diff.mean, phi.diff.q2.5, phi.diff.q97.5, phi.beta.mean, phi.beta.q2.5, phi.beta.q97.5, deparse.level = 1)

write.csv(Simulation.results,"Simulation_results.csv")

Simulation.all <- cbind(Simulation,Simulation.results)
write.csv(Simulation.all,"Simulation_params and results.csv")





