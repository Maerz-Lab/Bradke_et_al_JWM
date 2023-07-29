# Bradke. D. R., B. A. Crawford, M. Kaylor, and J. C. Maerz. 2023. Facing uncertainty: evaluating and improving a standard monitoring method for diamond-backed terrapins. Journal of Wildlife Management.

# Code to fit an exponential decay function to percent of simulations with survival increase detected vs number of sites sampled for each certainty level. 

sites <- rep(c(1,2,4,6,8),4) # numbers of sites considered
confidence <- rep(c(95,85,75,"mean"), each=5) #levels of certainty considered 
detected <- c(34,50,79,92,96,61,74,91,95,99,75,80,96,99,100,90,96,100,100,100) #percent of simulations with survival increase detected

df <- as.data.frame(cbind(sites, confidence, detected)) #create data frame

#Format data for the model
df$detected <- as.numeric(df$detected)
df$sites <- as.numeric(df$sites)
df.wide <- reshape(df, idvar = "sites", timevar = "confidence", direction = "wide")
y <- as.matrix(df.wide[,-1])
site <- df.wide$sites

####################
#### Model code ####
####################

library(jagsUI)

cat(file = "mod_jags.txt","
model {

#priors

for (c in 1:4) {
sd[c] ~  dunif(0,10)
tau[c] <- 1 / (sd[c] * sd[c]) #precision

y0[c] ~ dunif(0,100) 
log.alpha[c] ~ dunif(-1.5,0) #peak time
}

for (c in 1:4) {
for (i in 1:5) {
y[i,c] ~ dnorm(mean.y[i,c],tau[c]) #y = the data supplied to model

mean.y[i,c] <- 100 + (y0[c]-100)*exp(-(exp(log.alpha[c])*(site[i]-1))) #exponential decay function 

}
}

}
")

# bundle data
bdata <- list(y=y,site=site)

# initial values
inits <- function() list(y0 = runif(n=4,min=0,max=25), log.alpha = runif(n=4,min=-1,max=0), sd = runif(n=4,min=0,max=3))

# params to follow
params <- c("y0", "log.alpha", "sd", "mean.y")

# MCMC settings
na <-20000  ; ni <- 120000 ; nt <- 1 ; nb <- 10000 ; nc <- 3

# run model
out <- jags(bdata, inits, params, "mod_jags.txt", n.adapt = na,
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# view results
print(out, 3)


