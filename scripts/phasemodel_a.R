# Load Library and Data ----
library(here)
library(dplyr)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)

`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','c14data_farmers.RData'))

# General Setup ----
# Data --
dat <- list(cra = f_dateInfo$cra,
            cra_error = f_dateInfo$cra_error,
            constraint_uniform = rep(1, f_constants$n_areas),
            cra_constraint = rep(1, f_constants$n_dates)) # Set-up constraint for ignoring inference outside calibration range

# Initial parameters --
buffer <- 100
theta_init <- f_dateInfo$median_dates
delta_init <- f_siteInfo$diff + buffer
alpha_init <- f_siteInfo$earliest + buffer/2


#Calibration curve
f_constants$cc <- as.numeric(as.factor(f_dateInfo$calCurve)) #intcal20==1 and shcal20==2

# Dummy extension of the calibration curve
f_constants$calBP <- c(1000000, f_constants$calBP, -1000000)
f_constants$C14BP <- rbind(c(1000000,1000000), f_constants$C14BP, c(-1000000,-1000000))
f_constants$C14err <- rbind(c(1000,1000), f_constants$C14err, c(1000,1000))

# Initialise regional parameters ----
# Initialise sq areas which contain sites
init_a  <- aggregate(earliest~area_id, FUN=max, data=f_siteInfo) #find earliest date in each region k
init_b  <- aggregate(latest~area_id, FUN=min, data=f_siteInfo) #find latest date in each area k

#Initialise sq areas which do not contain sites
init_empty_area <- function(init_df) {
  for(i in 1:f_constants$n_area){
    
    area_ids <- init_df$area_id #List of sq areas ids with sites
    
    if (i %!in% area_ids){
      empty_sq_id <- i #Id of empty hex
      neighbour_sq <- which.min(abs(i - area_ids)) #Determine closest sq neighbor which has sites. If there are more than one neighbour sq area with sites, it selects the first observation (i.e. sq with the smallest id, since ids are in ascending order -- which favours neighbours South West to the sq area in question (proving a good lower bound to age estimate). This makes sense, since our putative origin is in Sq 12 and we think the direction of spread is North East). #TODO: refine this selection of neighbours, once general direction of spread is determined through GPQR 
      neighbour_sq_id <- area_ids[neighbour_sq] #Determine area id of closest sq neighbor
      neighbour_date <- init_df[neighbour_sq , 2] #Select the date associated with the neighbour sq 
      init_df <- rbind(init_df, c(i, neighbour_date)) #Assign this date to the empty sq area
    }
  }
  return(init_df)
}
init_a <- init_a %>% init_empty_area() %>%  arrange(area_id)
init_b <- init_b %>% init_empty_area() %>%  arrange(area_id)

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

#===============================================================================
# MCMC RunScript (Uniform Model a) ----
unif_model_a  <- function(seed, d, theta_init, alpha_init, delta_init, init_a, init_b, constants, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
      alpha[j] ~ dunif(max = a[id_area[j]], min = b[id_area[j]]);
    }

    for (i in 1:n_dates){
      theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
      #Calibration
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ , cc[i]]); #c14age #Index cc selects the correct calibration curve
      cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ , cc[i]]);
      sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=mu[i], sd=sd[i]);
    }

    # Set Prior for Each Region
    for (k in 1:n_areas){
      a[k] ~ dunif(4000, 12000);
      b[k] ~ dunif(4000, 12000);
      constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    }
    # Hyperprior for duration
    gamma1 ~ dunif(1,20) #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500) #Hyperprior for mode
  })
  
  # Define Initial values ----
  inits <- list(theta=theta_init, 
                alpha=alpha_init, 
                delta=delta_init, 
                a=init_a, 
                b=init_b)
  inits$gamma1  <- 10
  inits$gamma2  <- 200
  
  # Compile and Run model	----
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf$addMonitors(c('theta','delta','alpha'))
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) 
}

# Run MCMCs ----

# MCMC Setup
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12, 34, 56, 78)
niter  <- 6000000
nburnin  <- 3000000
thin  <- 300

out_unif_model_a  <-  parLapply(cl = cl, 
                                X = seeds, 
                                fun = unif_model_a, 
                                d = dat, 
                                constants = f_constants, 
                                theta_init = theta_init, 
                                alpha_init = alpha_init, 
                                delta_init = delta_init,  
                                init_a = init_a, 
                                init_b = init_b, 
                                niter = niter, 
                                nburnin = nburnin,
                                thin = thin)

out_unif_model_a <- mcmc.list(out_unif_model_a)


# Diagnostics ----
rhat_unif_model_a <- gelman.diag(out_unif_model_a, multivariate = FALSE)
ess_unif_model_a <- effectiveSize(out_unif_model_a)
agg_unif_model_a <- agreementIndex(dat$cra,
                                   dat$cra_error,
                                   calCurve = f_dateInfo$calCurve,
                                   theta = out_unif_model_a[[1]][ , grep("theta", colnames(out_unif_model_a[[1]]))],
                                   verbose = F)


#-------------------------------------------------------------------------------
# Save output ----
save(out_unif_model_a, 
     rhat_unif_model_a, 
     ess_unif_model_a, 
     agg_unif_model_a, 
     file=here("output","phase_model_a.RData"))
