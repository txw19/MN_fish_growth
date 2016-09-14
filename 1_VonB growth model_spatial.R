# rm(list=ls())

#### Load libraries
library(MCMCpack) # rwish function
library(R2jags)
require(lubridate)

# Read in data
dat<-read.csv('MN_age_length_WAE_YEP_LMB.csv',na.strings='NA', header=T)
head(dat)
dim(dat)

unique(dat$YEAR)
unique(dat$SPECIES_NAME)

# Select Species
dat <- dat[dat$SPECIES_NAME=='Walleye',]
dim(dat)

# Get rid of unused level
dat$SPECIES_NAME <- droplevels(dat$SPECIES_NAME)

# Grab useful columns
dat <- dat[,c('FISHERIES_WATERBODY_ID','ASSIGNED_AGE_YEARS','TOTAL_LENGTH_MM','SURVEY_ID_DATE')]
head(dat)

dat$date <- strptime(dat$SURVEY_ID_DATE, "%e-%b-%y") 

head(dat)

#Note:  These extractions require the lubridate package to be installed
# Extract year
dat$year <- year(dat$date)
head(dat)

# Get rid of 200 year old fish
dat <- dat[dat$ASSIGNED_AGE_YEARS<100,]

# Keep only recentyears of data
dat <- dat[dat$year >=2010,]

dim(dat)


###########################################################
## Obtain number of observations per site, called length.n
#Activate downloaded library
library(doBy)
sumfun <- function(x, ...){
  c(n=length(x))
}
#Run summary of data
sum1<-summaryBy(TOTAL_LENGTH_MM~FISHERIES_WATERBODY_ID,data=dat,FUN=sumfun,na.rm=T) 
dim(sum1)
head(sum1)
min(sum1$TOTAL_LENGTH_MM.n)
max(sum1$TOTAL_LENGTH_MM.n)
head(sum1)

####### Keep lakes with > N observations
sum2 <- sum1[which(sum1$TOTAL_LENGTH_MM.n > 75), ]
dim(sum2)
head(sum2)

dat <- dat[dat$FISHERIES_WATERBODY_ID%in%sum2$FISHERIES_WATERBODY_ID,]
dim(dat)

# Number of lakes
length(unique(dat$FISHERIES_WATERBODY_ID))


# Quick and ugly plots
plot(TOTAL_LENGTH_MM ~ jitter(ASSIGNED_AGE_YEARS), data=dat,pch=16,xlab='Age (yrs)',ylab='Length (mm)')

head(dat)

# sort by waterbody_id
dat <- dat[order(dat$FISHERIES_WATERBODY_ID),]
head(dat)


#################################################################
########## BUGS CODE ############################################
#################################################################


sink("vonBmodel.txt")
cat("
    model{
    for(i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- Linf[g[i]] * (1-exp(-k[g[i]] * (age[i] - t0[g[i]] )))
    }
    
    tau.y <- pow(sigma.y,-2)
    sigma.y ~ dunif(0,100)

# Parameters modeled on log-scale
for(j in 1:J){
  Linf[j] <- exp(BB[j,1])
	k[j] <- exp(BB[j,2])
	t0[j] <- exp(BB[j,3])-10 # A constant of 10 is added (subtracted) to t0 to ensure that negative values are possible, becuase t0 is estimated on log-scale
	BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,]) # Multivariate normal dist'n;  Tau.B is a precision 
	BB.hat[j,1] <- mu.Linf
	BB.hat[j,2] <- mu.k
	BB.hat[j,3] <- mu.t0
}

  # Priors for population-average parameters
  mu.Linf ~ dnorm(0,.0001)
  mu.k ~ dnorm(0,.0001)
  mu.t0 ~ dnorm(0,.0001)


# Model variance-covariance
  Tau.B[1:K,1:K] ~ dwish(W[,], df)
  df <- K+1
  Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
  for (k in 1:K){
    for (k.prime in 1:K){
      rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
        sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
  }

    
} # end model
    ",fill=TRUE)
sink()

################
# Number of sites
J <-length(unique(dat$FISHERIES_WATERBODY_ID))

# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 3

# Create identity matrix for Wishart dist'n
W <- diag(3)

dat$g <- as.factor(as.numeric(as.factor(dat$FISHERIES_WATERBODY_ID)))

# load data
data <- list(y = dat$TOTAL_LENGTH_MM, age = dat$ASSIGNED_AGE_YEARS, g = dat$g, n = dim(dat)[1],
             J = J, W=W, K=K )


# Initial values
inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),
                         sigma.y = runif(1,1,10), 
                    BB=array(c(rep(log(600) +rnorm(1,0.01,0.01),J),rep(log(0.4)+rnorm(1,0.001,0.1),J),rep(log(0.5+10)+rnorm(1,0.01,0.1),J)),
                                 c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }



# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B")


# MCMC settings
ni <- 50000
nt <- 3
nb <- 40000
nc <- 3


############################################
###### DO analysis in JAGS on single core
start.time = Sys.time()         # Set timer (7.67 mins)

out <- jags(data = data, inits = inits, parameters.to.save = params1, 
            model.file = "vonBmodel.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize the result
print(out, digits = 2)
# str(out)

# Find which parameters, if any, have Rhat > 1.1
which(out$BUGSoutput$summary[, c("Rhat")] > 1.1)

