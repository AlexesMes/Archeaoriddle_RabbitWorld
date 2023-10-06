# Load Libraries and Data ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(parallel)
library(microbenchmark)

#-------------------------------------------------------------------------------
## Data Setup ----
# Load and prepare data ----
load(here('data','c14data_farmers.RData'))

## Define Constants and Data ----
dateInfo  <- subset(f_dateInfo, earliestAtSite==TRUE)
dateInfo  <- dateInfo[order(dateInfo$siteID, decreasing = FALSE),] #Arrange in ascending site ID order

# Generate list of observed data
dat  <- list(cra=dateInfo$cra, 
             cra_error=dateInfo$cra_error)

# Constants
f_constants$tau <- 0.90

#Calibration curve
f_constants$cc <- as.numeric(as.factor(dateInfo$calCurve)) #intcal20==1 and shcal20==2

# Dummy extension of the calibration curve
f_constants$calBP <- c(1000000, f_constants$calBP, -1000000)
f_constants$C14BP <- rbind(c(1000000,1000000), f_constants$C14BP, c(-1000000,-1000000))
f_constants$C14err <- rbind(c(1000,1000), f_constants$C14err, c(1000,1000))

# Fixed Inits
theta_init <- dateInfo$median_dates

#-------------------------------------------------------------------------------
## Define Run Function ----
runFun  <- function(seed,dat, theta_init, constants, niter, nburnin, thin)
{
  library(nimbleCarbon)
  library(truncnorm)
  library(cascsim)
  
  # Variance Covariance Function ----
  cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0), sigmasq = double(0)) 
  {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE) #Results matrix
    deltaij <- matrix(nrow = n, ncol = n, init = TRUE) #Identity matrix
    diag(deltaij) <- 1
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
    return(result)
  })
  Ccov_ExpQ <- compileNimble(cov_ExpQ)
  assign('cov_ExpQ', cov_ExpQ, envir=.GlobalEnv)
  
  # Handle constraints on dispersal rate
  dat$lim  <- rep(1, constants$n_sites)
  # Handle constraints on calibration so that everything is within calibration curve
  dat$cra_constraint = rep(1, constants$n_dates)
  
  # Core model ----
  model <- nimbleCode({
    for (i in 1:n_sites){
      # Model
      rate[i] <- -1/(s[i]-beta1) #s_i allows for site-level variation in the dispersal rate
      lim[i] ~ dconstraint(rate[i]>0)
      mu[i] <- beta0 + (s[i]-beta1)*dist_org[i]
      theta[i] ~ dAsymLaplace(mu=mu[i], sigma=sigma, tau=tau)
      #Calibration
      mu.date[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ , cc[i]]); #c14age #Index cc selects the correct calibration curve
      cra.constraint[i] ~ dconstraint(mu.date[i] < 50193 & mu.date[i] > 95)
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ , cc[i]]);
      sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=mu.date[i], sd=sd[i]);
    }
    #priors
    beta0 ~ dnorm(8000, sd=200);
    beta1 ~ dexp(1)
    sigma ~ dexp(0.01)
    etasq ~ dexp(20);
    rho ~ T(dgamma(10, (10-1)/100), min=1, max=620); #Rho is the length scale parameter defining the rate the covariance declines as a function of distance, i.e. it controls the extent of the spatial autocorrelation, currently defined to be uniform, with a mode 150 and defined between 1km and 620km (the smallest and largest inter-distance between sampled sites). #TODO: build in a varying spatial autocorrelation -- see notes.
    mu_s[1:n_sites] <- 0;
    cov_s[1:n_sites, 1:n_sites] <- cov_ExpQ(dist_mat[1:n_sites, 1:n_sites], rho, etasq, sigmasq=0.000001) #Covariance Matrix. sigmasq is assigned a small positive constant value (10^-6) for programmatic reasons (non-positive eigenvectors). This parameter could, however, be used to introduce additional intra-site covariance, if necessary. 
    s[1:n_sites] ~ dmnorm(mu_s[1:n_sites], cov = cov_s[1:n_sites, 1:n_sites])
  }) 
  
  # MCMC initialisation ----
  set.seed(seed)
  inits  <-  list()
  inits$theta  <- theta_init
  inits$beta0 <- rnorm(1, 8000, 200)
  inits$beta1 <- rexp(1, rate=2)
  inits$sigma  <- rexp(1, 0.01)
  inits$rho  <- rtgamma(1, shape=10, scale=(10-1)/100, min=1, max=620)
  inits$etasq  <- rexp(1,20)
  inits$s  <- rep(0, constants$n_sites) #create vector
  inits$cov_s <- Ccov_ExpQ(constants$dist_mat, inits$rho, inits$etasq, 0.000001)
  inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$n_sites) #chol() performs the Cholesky decomposition, t() transpose
  inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix
  
  # Model Compilation ----
  model.gpqr <- nimbleModel(model, constants = constants, data=dat, inits=inits)
  cModel.gpqr <- compileNimble(model.gpqr)
  
  # MCMC configuration ----
  conf.gpqr <- configureMCMC(model.gpqr, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf.gpqr$addMonitors('s')
  conf.gpqr$addMonitors('rho')
  conf.gpqr$addMonitors('etasq')
  conf.gpqr$removeSamplers('s[1:43]') #select relevant vector, where n=43 is the number of sites 
  conf.gpqr$removeSamplers('beta1')
  conf.gpqr$addSampler(c('beta1','s[1:43]'), type='AF_slice') 
  conf.gpqr$removeSamplers('beta0');conf.gpqr$addSampler('beta0', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  conf.gpqr$removeSamplers('sigma');conf.gpqr$addSampler('sigma', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  conf.gpqr$removeSamplers('etasq');conf.gpqr$addSampler('etasq', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  conf.gpqr$removeSamplers('rho');conf.gpqr$addSampler('rho', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  MCMC.gpqr <- buildMCMC(conf.gpqr)
  cMCMC.gpqr <- compileNimble(MCMC.gpqr)
  
  # MCMC execution
  results <- runMCMC(cMCMC.gpqr, nchain=1, niter = niter, thin=thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) 
  return(results)
}


## Setup and Execution of MCMC in Parallel ----
ncores <- 4
cl <- makeCluster(ncores)

# Run the model in parallel:
seeds <- c(12, 45, 67, 89)
niter = 6000000
nburnin = 3000000
thin = 300

#Run the model in parallel
chain_output <- parLapply(cl = cl, X = seeds, fun = runFun, dat = dat, constants = f_constants, theta = theta_init, niter = niter, nburnin = nburnin, thin = thin)
stopCluster(cl)

# Convert into a mcmc.list object for diagnostic (see below)
gpqr_tau90 <- coda::mcmc.list(chain_output)
rhat.t90 <- coda::gelman.diag(gpqr_tau90, multivariate = F)
# range(rhat.t90$psrf[,1])
# ii  <- which(rhat.t90$psrf[,1]>1.01)  
ess.t90  <- coda::effectiveSize(gpqr_tau90)
# range(ess.t90)


#-------------------------------------------------------------------------------
## Store Output ----
save(rhat.t90, ess.t90, gpqr_tau90, file=here('output', 'gpqr_tau90.RData'))