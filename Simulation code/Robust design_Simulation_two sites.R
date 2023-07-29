# Bradke. D. R., B. A. Crawford, M. Kaylor, and J. C. Maerz. 2023. Facing uncertainty: evaluating and improving a standard monitoring method for diamond-backed terrapins. Journal of Wildlife Management.

# The "SIMULATION CODE" is for simulating diamond-backed terrapin capture-recapture datasets with a positive change in survival. 
# The "MODEL CODE" is used to evaluate if the simulated increases in survival would be detected using standard monitoring methods if data were collected at two sites.

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

CH <- lapply(1:100, matrix, data= 0, nrow=NT, ncol=N2) 

tide.site <- array(dim=c(iter,N1-1)) 

#to draw random MCMC iterations and associated parameter values from the RD model results to be used in simulation iterations
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

write.csv(Simulation,"Simulation_params2.csv")

tide.actual <- cov$Actual_extrap_Jekyll #bring in actual tide vals from sampling at Jekyll during years 1-10 and 12-13 (no sampling in year 11)
tide.mean <- rep(NA,12)
cnss.actual <- c(0,3,6,9,12,15,18,21,24,27,30,32,34)
for(i in 1:10){
  tide.mean[i] <- mean(tide.actual[(cnss.actual[i]+1):(cnss.actual[i]+3)])
}
for(i in 12:13){
  tide.mean[i-1] <- mean(tide.actual[(cnss.actual[i]+1):(cnss.actual[i]+2)])
}


for (k in 1:iter){
  logit.phi <- c(rnorm(pre.int,qlogis(pre.phi),phi.sd[k]),rnorm(post.int,qlogis(post.phi),phi.sd[k]))
  phi <- plogis(logit.phi) #phi for each sampling interval
  
  logit.p <- rnorm(N2,qlogis(p.mean[k]),p.sd[k])
  
  tide.site[k,] <- rnorm((N1-1), mean(tide.mean),sd(tide.mean)) #select tide covariate vals for simulation
  
  logit.gam <- qlogis(gam.mean[k])+gam.beta[k]*tide.site[k,]
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
    
    exist.occ <- c(rep(1,NS), rep(2:N1,each = NN)) #vector of occasions when individuals first "exist"
    
    # Fill the CH matrix
    for (i in 1:NT){
      
      #simulate if captured in first secondary occasion 
      #of first primary period that individual exists
      rp <- rbinom(1, 1, plogis(LOGIT.P[i,(cnss[exist.occ[i]]+1)]))
      if (rp==1) {CH[[k]][i,(cnss[exist.occ[i]]+1)] <- 1}
      #simulate if captured in remaining secondary occasions 
      #of first primary period that individual exists
      for (j in 2:nss[exist.occ[i]]){
        P <- plogis(LOGIT.P[i,(cnss[exist.occ[i]]+j)]+C.BETA*max(CH[[k]][i,(cnss[exist.occ[i]]+1):(cnss[exist.occ[i]]+j-1)]))
        rp <- rbinom(1, 1, P)
        if (rp==1) {CH[[k]][i,(cnss[exist.occ[i]]+j)] <- 1}
      }#j
      
      if (exist.occ[i]==N1) {next} #If first exists in final primary period, skip to next individual
      for (t in (exist.occ[i]+1):N1){
        sur <- rbinom(1, 1, PHI[i,(t-1)]) #did individual survive to second primary period?
        if (sur==0) {break}	# If dead, move to next individual 	
        fid <- rbinom(1, 1, GAM[i,(t-1)]) #is individual available in second primary period?
        if (fid==0) {next}
        #if available, was individual captured in first secondary occasion?
        rp <- rbinom(1, 1, plogis(LOGIT.P[i,(cnss[t]+1)]))
        if (rp==1) {CH[[k]][i,(cnss[t]+1)] <- 1}         
        for (j in 2:nss[t]){
          #Was individual captured in eaCH[[k]] remaining secondary occasion? 
          P <- plogis(LOGIT.P[i,(cnss[t]+j)]+C.BETA*max(CH[[k]][i,(cnss[t]+1):(cnss[t]+j-1)]))
          rp <- rbinom(1, 1, P)
          if (rp==1) {CH[[k]][i,(cnss[t]+j)] <- 1}
        }#j
      } #t
    } #i
    return(CH[[k]])
  }
  
  # Execute function
  CH[[k]] <- simul.rd(PHI, LOGIT.P, C.BETA, GAM, NT, nss, cnss) 
  
  CH[[k]] <- CH[[k]][rowSums(CH[[k]])>0, ] #delete individuals with all "0" capture histories
  
}


n.sites <- 2
iter_sim <- 100
start_sim <- 1
end_sim <- 100

ch.sites.2 <- array(dim=c(iter*n.sites))
ch.sites.2 <- as.data.frame(ch.sites.2)
cnss2 <- seq(0,iter*n.sites,by=n.sites)

CH.comb <- vector(mode='list', length=100)
tide <- vector(mode='list', length=100)
group <- vector(mode='list', length=100)

for (k in 1:iter){ 
  ch.sites <- sample(1:100,n.sites,replace=F)
  ch.sites.2[(cnss2[k]+1):cnss2[k+1],] <- ch.sites
  ch.sites.2[(cnss2[k]+1),2:16] <- Simulation[ch.sites[1],]
  ch.sites.2[cnss2[k+1],2:16] <- Simulation[ch.sites[2],]
  
  CH1 <- as.data.frame(CH[[ch.sites[1]]])
  CH2 <- as.data.frame(CH[[ch.sites[2]]])
  
  group[[k]] <- c(rep(1,dim(CH1)[1]),rep(2,dim(CH2)[1]))
  
  CH.comb[[k]] <- as.matrix(rbind(CH1,CH2))
  
  tide[[k]] <- rbind(tide.site[ch.sites[1],],tide.site[ch.sites[2],])
}


#create file names to store simulation results later
file.num <- c(start_sim:end_sim)

file.samps <- paste("./Samples2/Simulation_samples",file.num,".csv",sep = "")
file.summs <- paste("./Summary2/Simulation_summary",file.num,".csv",sep = "")

phi.diff.mean <- array(dim=c(iter_sim))
phi.diff.q2.5 <- array(dim=c(iter_sim))
phi.diff.q97.5 <- array(dim=c(iter_sim))


  #############################################
  ################ MODEL CODE##################
  #############################################

    
for (k in start_sim:end_sim){     
  #To get matrix of capture-recapture data for primary periods only
  z.data <- matrix(0, ncol = N1, nrow = dim(CH.comb[[k]])[1])
  for (i in 1:dim(CH.comb[[k]])[1]){
    for (t in 1:N1){
      cap <- sum(CH.comb[[k]][i,(cnss[t]+1):(cnss[t]+per.prime)])
      if (cap>0) z.data[i,t] <- 1
    } #t
  } #i
  
  #To get "if captured before" per primary period
  X <- array(0,dim=dim(CH.comb[[k]]))
  for (i in 1:dim(CH.comb[[k]])[1]){
    for (t in 1:N1){
      for (j in 2:nss[t]) { 
        if(CH.comb[[k]][i,(cnss[t]+j-1)]==1) X[i,(cnss[t]+j)] <- 1
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
  
  BRD <- c(rep(1,pre.int),rep(2,post.int)) #Create binary temporal co-variate indicating pre (0) vs post(1) regulations
  
  
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
        logit(phi[i,t]) <- phi.mu[BRD[t]] + phi.eps[group[i],t] 
      } #t
    } #i
    
    for (u in 1:g){
    for (t in 1:(n1-1)){
      phi.eps[u,t] ~ dnorm(0,phi.tau[u])
    } #t
    }
    
    for (b in 1:2){
    phi.mean[b]  ~ dunif(0, 1)       #Priors for mean phi - 1 for pre and one for post BRD
    phi.mu[b] <- log(phi.mean[b]/(1-phi.mean[b]))     
    }

    for (u in 1:g){
    phi.tau[u] <- pow(phi.sd[u],-2)
    phi.sd[u] ~ dunif(0,5) #Prior for sd (random year effect) on survival 
    }
    
    #Availability
    
    for (i in 1:n) {
      for (t in f[i]:(n1-1)){
        logit(gamma[i,t]) <- gam.mu[group[i]] + gam.beta[group[i]]*tide.mean[group[i],t]
      } #t
    } #i
    
    for (u in 1:g){
    gam.mean[u] ~ dunif(0,1) #Prior for mean availability
    gam.mu[u] <- log(gam.mean[u]/(1-gam.mean[u])) 
    gam.beta[u] ~ dnorm(0,0.37) #Prior for fixed tide effect on availability
    }
    
    #Capture parameters
    for (u in 1:g){
    for (j in 1:n2){
      p.eps[u,j] ~ dnorm(0,p.tau[u])
    } #j
    }
    
    for (i in 1:n) {
      for (j in 1:n2) {
        logit(p[i,j]) <- p.mu[group[i]] + p.eps[group[i],j] + c.beta[group[i]]* X[i,j] 
      } #j
    } #i
    
    for (u in 1:g){
    p.mean[u] ~ dunif(0,1) #Prior for mean detection
    p.mu[u] <- log(p.mean[u]/(1-p.mean[u])) 
    p.tau[u] <- pow(p.sd[u],-2)
    p.sd[u] ~ dunif(0,5) #Prior for sd (random day effect) on capture probability 
    c.beta[u] ~ dnorm(0, 0.37) #Prior for within secondary period recapture effect (i.e., trap shy effect)
    }
    
    
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
    
    phi.diff <- phi.mean[2]-phi.mean[1] #effect size of BRD reg on survival
    
  }) #End Model
  
  #constants
  Mod6Consts <- list(f = f, n=nrow(CH.comb[[k]]),n1=N1,n2=N2,nss=nss,cnss=cnss, g =length(unique(group[[k]])))
  #data
  Mod6Data <- list(y=CH.comb[[k]],tide.mean=tide[[k]],z=z.known, BRD=BRD, s=s.known,X=X, group=group[[k]])
  
  #Parameters monitored
  params<- c("phi.mean","phi.sd","phi.diff","p.mean","p.sd","c.beta","gam.beta","gam.mean")
  
  #Set initial values
  inits <- list(phi.mean=runif(length(unique(BRD)),0.5,1),phi.sd=runif(length(unique(group[[k]])),0,5),phi.beta=rnorm(1),gam.mean=runif(length(unique(group[[k]])),0,1),gam.beta=rnorm(length(unique(group[[k]])),0,1),p.mean=runif(length(unique(group[[k]])),0,1),c.beta=rnorm(length(unique(group[[k]])),0,1),p.sd=runif(length(unique(group[[k]])),0,5),z=z.inits)
  
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
    nimbleOptions(showCompilerOutput = TRUE)
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
  
  phi.diff.mean[k] <- summary['phi.diff',1]
  phi.diff.q2.5[k] <- summary['phi.diff',3]
  phi.diff.q97.5[k] <- summary['phi.diff',5]
  
}

Simulation.results <- cbind(phi.diff.mean, phi.diff.q2.5, phi.diff.q97.5, deparse.level = 1)

write.csv(Simulation.results,"Simulation_results2.csv")

colnames(ch.sites.2) <- c("sim",colnames(Simulation))
write.csv(ch.sites.2,"Simulation_CH_params2.csv")








