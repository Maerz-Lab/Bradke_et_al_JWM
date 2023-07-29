# Bradke. D. R., B. A. Crawford, M. Kaylor, and J. C. Maerz. 2023. Facing uncertainty: evaluating and improving a standard monitoring method for diamond-backed terrapins. Journal of Wildlife Management.

# This code is for running the set of robust-design CSJ candidate models on the capture-recapture data collected from two salt marsh tidal creeks adjacent to the Downing Musgrove Causeway in Glynn County, Georgia during 2010 – 2022 (except for 2020) 
# Individual models are described within the model text files

# load required packages
library(jagsUI)
library(coda)
library(loo) 
library(beepr)
library(MCMCvis)

#Read in data files
terpRD <- read.csv("RD_Data_2010_2022_Males.csv") # capture history data from seining (2010-2022)
cov <- read.csv("Tide_cov_2010_2022.csv") # Tide covariate 
z.data <- read.csv("Primary Per Data for z.csv") # Primary period capture histories to produce known z states 
DS <- read.csv("Dummy_Secondary.csv") # Binary covariate indicating if survey was conducted (1) or not (0) to account for no surveys in 2020 (for capture probability)

#Sampling information and covariates
n1=13 # number of primary periods
n2=36 # number of total secondary sampling periods
nss=c(3,3,3,3,3,3,3,3,3,3,2,2,2) # number of secondary periods in primary period i
cnss=c(0,3,6,9,12,15,18,21,24,27,30,32,34) # cumulative sum of nss from 1 to i-1
# extract tide value for each secondary occasion - daily covariate for capture probability
tide=cov$Actual_extrap_Jekyll 
# calculate mean tide within each primary period - annual covariate for availability
tide.mean <- rep(NA,n1)
for(i in 1:10){
  tide.mean[i] <- mean(tide[(cnss[i]+1):(cnss[i]+3)])
}
for(i in 12:13){
  tide.mean[i] <- mean(tide[(cnss[i]+1):(cnss[i]+2)])
}
tide.mean[11] <- mean(tide.mean,na.rm = T) # use mean of other values because no sampling this year

tide.mean <- tide.mean[-1] # remove first covariate value because availability not estimated for year 1

#remove column with IDs
terpRD <- as.matrix(terpRD[,2:37]) 
z.data <-as.matrix(z.data[,2:14])
DS <- as.matrix(DS[,2:37])

#To get "if captured before" per primary period - covariate for recapture effect
X <- array(0,dim=dim(terpRD))
for (i in 1:dim(terpRD)[1]){
  for (t in 1:n1){
    for (j in 2:nss[t]) { 
      if(terpRD[i,(cnss[t]+j-1)]==1) X[i,(cnss[t]+j)] <- 1
      if(X[i,(cnss[t]+j-1)]==1) X[i,(cnss[t]+j)] <- 1
    }
  }
}

# To get occasion of first capture for each individual
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

#function for initial values for z state
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}
# initial values for z state
z.inits <- cjs.init.z(z.data,f)

# known states for availability based on capture histories
s.known <- as.matrix(z.data)
s.known[s.known==0]<-NA

#MCMC settings
ni<-120000
nt<-1
nb<-30000
nc<-3


####Run Mod1

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,tide=tide,z=z.known,s=s.known, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean", "alpha", "gam.mean", "gam.sd","loglike.new", "phi.t", "p.t")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),gam.sd=runif(1,0,5),p.mean=runif(1,0,1),alpha=rnorm(1),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod1 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod1.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod1)

pdf(file="Traceplots_mod1.pdf", height=11, width=8.5)  # for boxier images
MCMCtrace(jmod1,
          params = c("phi.mean", "phi.sd","p.mean", "alpha", "gam.mean", "gam.sd"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod1 <- jmod1$summary
write.csv(summary.mod1, "Mod1_summary.csv")

joint_Loglik_mod1 <- jmod1$sims.list$loglike.new
joint_waic_m1 <- waic(joint_Loglik_mod1)
joint_waic_m1$estimates["waic",]


####Run Mod2

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,tide=tide,z=z.known,s=s.known,X=X, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean", "delta","alpha", "gam.mean", "gam.sd","loglike.new","phi.t", "time.p", "time.c")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),gam.sd=runif(1,0,5),p.mean=runif(1,0,1),delta=rnorm(1),alpha=rnorm(1),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod2 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod2.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod2)
pdf(file="Traceplots_mod2.pdf", height=11, width=8.5)  # for boxier images 
MCMCtrace(jmod2, 
          params = c("phi.mean", "phi.sd","p.mean", "alpha", "gam.mean", "gam.sd"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod2 <- jmod2$summary
write.csv(summary.mod2, "Mod2_summary.csv")

joint_Loglik_mod2 <- jmod2$sims.list$loglike.new
joint_waic_m2 <- waic(joint_Loglik_mod2)
joint_waic_m2$estimates["waic",]


####Run Mod3

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,tide.mean=tide.mean, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean", "p.sd","gam.mean", "beta","loglike.new", "phi.t", "p.t")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),beta=rnorm(1),p.mean=runif(1,0,1),p.sd=runif(1,0,5),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod3 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod3.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod3)
pdf(file="Traceplots_mod3.pdf", height=11, width=8.5)  # for boxier images 
MCMCtrace(jmod3, 
          params = c("phi.mean", "phi.sd","p.mean", "p.sd","gam.mean", "beta"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod3 <- jmod3$summary
write.csv(summary.mod3, "Mod3_summary.csv")

joint_Loglik_mod3 <- jmod3$sims.list$loglike.new
joint_waic_m3 <- waic(joint_Loglik_mod3)
joint_waic_m3$estimates["waic",]


####Run Mod4

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,X=X,tide.mean=tide.mean, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean","p.sd","gam.mean","beta","loglike.new", "delta", "time.p", "time.c", "phi.t")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),beta=rnorm(1),delta=rnorm(1),p.mean=runif(1,0,1),p.sd=runif(1,0,5),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod4 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod4.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod4)

pdf(file="Traceplots_mod4.pdf", height=11, width=8.5)  # for boxier images 
MCMCtrace(jmod4, 
          params = c("phi.mean", "phi.sd","p.mean", "p.sd","gam.mean", "beta","delta"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod4  <- jmod4$summary
write.csv(summary.mod4, "Mod4_summary.csv")

joint_Loglik_mod4 <- jmod4$sims.list$loglike.new
joint_waic_m4 <- waic(joint_Loglik_mod4)
joint_waic_m4$estimates["waic",]


####Run Mod5

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,tide.mean=tide.mean, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean", "p.sd","gam.mean", "beta","loglike.new","phi.t","p.t")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),beta=rnorm(1),p.mean=runif(1,0,1),p.sd=runif(1,0,5),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod5 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod5.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod5)

pdf(file="Traceplots_mod5.pdf", height=11, width=8.5)  # for boxier images 
MCMCtrace(jmod5, 
          params = c("phi.mean", "phi.sd","p.mean", "p.sd","gam.mean", "beta"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod5 <- jmod5$summary
write.csv(summary.mod5, "Mod5_summary.csv")

joint_Loglik_mod5 <- jmod5$sims.list$loglike.new
joint_waic_m5 <- waic(joint_Loglik_mod5)
joint_waic_m5$estimates["waic",]


####Run Mod6

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,X=X,tide.mean=tide.mean,Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean","p.sd","gam.mean","beta","loglike.new","delta")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),beta=rnorm(1),delta=rnorm(1),p.mean=runif(1,0,1),p.sd=runif(1,0,5),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod6 <- jags(data = rd.data, inits = inits, parameters.to.save = c(params, "phi.t", "time.p", "time.c", "pred.gamtide"), model.file = "mod6.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod6)

pdf(file="Traceplots_mod6.pdf", height=11, width=8.5)  # for boxier images 

MCMCtrace(jmod6, 
          params = c("phi.mean", "phi.sd","p.mean","p.sd","gam.mean","beta","delta"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()


summary_mod6  <- jmod6$summary
write.csv(summary_mod6, "Mod6_summary.csv")
samples <- jmod6$sims.list
saveRDS(samples, "Mod6.Rdata")
samples<- as.data.frame(samples)
write.csv(samples, "Mod6_samples.csv")

joint_Loglik_mod6 <- jmod6$sims.list$loglike.new
joint_waic_m6 <- waic(joint_Loglik_mod6)
joint_waic_m6$estimates["waic",]


####Run Mod6 - Last year of phi = fixed effect

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,X=X,tide.mean=tide.mean, Survey=DS)

#Parameters monitored
params <- c("phi.mean","phi.sd","p.mean","p.sd","gam.mean","alpha","beta","loglike.new","delta")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),alpha=rnorm(1),beta=rnorm(1),delta=rnorm(1),p.mean=runif(1,0,1),p.sd=runif(1,0,5),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod6_fixed <- jags(data = rd.data, inits = inits, parameters.to.save = c(params, "phi.t", "time.p", "time.c", "pred.gamtide"), model.file = "mod6_last year fixed.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod6_fixed)

pdf(file="Traceplots_mod6_fixed.pdf", height=11, width=8.5)  # for boxier images 

MCMCtrace(jmod6_fixed, 
          params = c("phi.mean", "phi.sd","p.mean","p.sd","gam.mean","beta","delta","alpha"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary_mod6_fixed <- jmod6_fixed$summary
write.csv(summary_mod6_fixed, "Mod6_fixed_summary.csv")
samples <- jmod6_fixed$sims.list
saveRDS(samples, "Mod6_fixed_samples.Rdata")
samples<- as.data.frame(samples)
write.csv(samples, "Mod6_fixed_samples.csv")

joint_Loglik_mod6_fixed <- jmod6_fixed$sims.list$loglike.new
joint_waic_m6_fixed <- waic(joint_Loglik_mod6_fixed)
joint_waic_m6_fixed$estimates["waic",]


####Run Mod7

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,tide.mean=tide.mean, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean", "gam.mean", "beta","loglike.new", "phi.t")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),beta=rnorm(1),p.mean=runif(1,0,1),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod7 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod7.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod7)

pdf(file="Traceplots_mod7.pdf", height=11, width=8.5)  # for boxier images 

MCMCtrace(jmod7, 
          params = c("phi.mean", "phi.sd","p.mean", "gam.mean", "beta"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod7  <- jmod7$summary
write.csv(summary.mod7, "Mod7_summary.csv")

joint_Loglik_mod7 <- jmod7$sims.list$loglike.new
joint_waic_m7 <- waic(joint_Loglik_mod7)
joint_waic_m7$estimates["waic",]


####Run Mod8

#Bundle data
rd.data<-list(y=terpRD, f = f, n=nrow(terpRD),n1=n1,n2=n2,nss=nss,cnss=cnss,z=z.known,s=s.known,X=X,tide.mean=tide.mean, Survey=DS)

#Parameters monitored
params <- c("phi.mean", "phi.sd","p.mean","delta","gam.mean","beta","loglike.new", "phi.t")

#Initial values
inits <- function() {list(phi.mean=runif(1,0.5,1),phi.sd=runif(1,0,5),gam.mean=runif(1,0,1),beta=rnorm(1),p.mean=runif(1,0,1),delta=rnorm(1),z=z.inits)}

#Run model
start.time <- Sys.time()
jmod8 <- jags(data = rd.data, inits = inits, parameters.to.save = params, model.file = "mod8.txt", n.chains=nc, n.adapt=NULL, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE, verbose=TRUE)
Sys.time() - start.time

beep("fanfare")

summary(jmod8)

pdf(file="Traceplots_mod8.pdf", height=11, width=8.5)  # for boxier images 

MCMCtrace(jmod8, 
          params = c("phi.mean", "phi.sd","p.mean","delta","gam.mean","beta"),
          ISB = FALSE, exact = TRUE, iter=90000,
          pdf = FALSE)
dev.off()

summary.mod8  <- jmod8$summary
write.csv(summary.mod8, "Mod8_summary.csv")

joint_Loglik_mod8 <- jmod8$sims.list$loglike.new
joint_waic_m8 <- waic(joint_Loglik_mod8)
joint_waic_m8$estimates["waic",]


############################################
# Produce WAIC table
############################################

WAIC_table <- loo_compare(joint_waic_m1, joint_waic_m2, joint_waic_m3, joint_waic_m4, joint_waic_m5, joint_waic_m6, joint_waic_m7, joint_waic_m8)[,1:8]
write.csv(WAIC_table,"WAIC_table.csv")



