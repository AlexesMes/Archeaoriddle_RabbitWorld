# Load Libraries and spatial data ----
library(here)
library(dplyr)
library(stringr)
library(truncnorm)
library(cascsim)
library(corrplot)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(maptools)
library(sf)
library(rgeos)
library(viridis)
library(cowplot)
library(wesanderson)
library(latex2exp)
library(gridExtra)
library(grid)
library(gridBase)
library(diagram)
library(quantreg)
library(coda)


`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data

# Load Observed Data
load(here('data','c14data_farmers.RData'))

# Load quantile regression results
load(here('output','quantreg_res.RData'))

#Read in basemap ----
rabbithole_map <- rast(here('maps', 'resources.tiff'))

#===============================================================================
# Generate Spatial Window for Analyses: Rabbithole ----
#Define sample window ----
boundary_win <- rabbithole_map %>% st_bbox() 
win.sf <- as.polygons(rabbithole_map > -Inf) %>% st_as_sf()
win <- as.polygons(rabbithole_map > -Inf) %>% as("Spatial")

#===============================================================================
##Bayesian Quantile Regression

## quantreg.R Figures ----

# Bayesian Quantile Regression Model with Measurement Error Plot ---- FIGURE 1
## Compute Fitted Model Confidence Intervals

# rq and median calibrated date
rq.ci <- predict.rq(fit_rq, newdata=data.frame(dist_org=0:400), interval='confidence') #Furthest site from proposed origin, Farwallow, is 364km

# Bayesian model 
qr.ch1 <- do.call(rbind,quantreg_sample)
post.alpha.quantreg <- qr.ch1[,'alpha']
post.beta.quantreg <- qr.ch1[,'beta']
post.alpha.beta.quantreg  <- data.frame(alpha=post.alpha.quantreg, beta=post.beta.quantreg)
post.theta.quantreg  <- qr.ch1[,grep('theta', colnames(qr.ch1))]
post.theta.med <- apply(post.theta.quantreg, 2, median)
post.quantreg <- apply(post.alpha.beta.quantreg, 1, function(x){x[1]-x[2]*0:400})
post.ci <- t(apply(post.quantreg, 1, quantile, c(0.025,0.5,0.975)))

## Plot and compare
# Transparency color utility function
col.alpha <- function(x,a=1){xx=col2rgb(x)/255;return(rgb(xx[1],xx[2],xx[2],a))}

pdf(file=here('output','figures','figure1.pdf'), width=8.5, height=7)
plot(NULL, xlim=c(0,400), ylim=c(8500,7000), axes=F, xlab='Distance from Farwallow Site (in km)', ylab='Cal BP') 
rect(xleft=-200, xright=4600, ybottom=2720, ytop=2350, col=col.alpha('grey',0.2), border=NA) #Demarcating Calibration Plateau Region
#axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)) #X-axis
#axis(2, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)) #Y-axis left
#axis(4, at=BCADtoBP(c(-1800, -1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=c('1800BC', '1400BC', '1000BC', '600BC', '200BC', '200AD', '600AD', '1000AD', '1400AD','1800AD'), cex.axis=0.6) #Y-axis right
points(f_constants$dist_org, f_siteInfo$earliest) #Median Calibrated Date
points(f_constants$dist_org, post.theta.med, pch=20) #Median Posterior theta
for (i in 1:nrow(f_siteInfo))
{
  lines(rep(f_constants$dist_org[i],2), c(f_siteInfo$earliest[i], post.theta.med[i]),lty=2)
}
lines(0:400, rq.ci[,1], lty=1, lwd=2, col='blue') #Quantile Regression on Median Dates
polygon(x=c(0:400, 400:0), c(rq.ci[,2], rev(rq.ci[,3])), col=col.alpha('lightblue', 0.4), border=NA)

lines(0:400, post.ci[,2], lty=1, lwd=2, col='indianred') #Bayesian Quantile Regression with Measurement Error
polygon(x=c(0:400, 400:0), c(post.ci[,1], rev(post.ci[,3])), col=col.alpha('indianred', 0.2), border=NA)


legend('bottomright', legend=c('Median Calibrated Date', TeX('Median Posterior $\\theta$'), 'Quantile Regression on Median Dates', 'Bayesian Quantile Regression with Measurement Error'), pch=c(1,20,NA,NA), lwd=c(NA,NA,2,2), col=c(1,1,'blue','indianred'), cex=0.8)
box()
dev.off()


#-------------------------------------------------------------------------------
# Posterior Dispersal Rate of non-spatial quantile regression Plot ---- FIGURE 2
pdf(here('output','figures','figure2.pdf'),height=5,width=5.5)
postHPDplot(1/post.beta.quantreg, xlab='km/year', ylab='Probability Density', prob=.90, main=TeX('Posterior of $1/\\beta_1$'))
dev.off()



#-------------------------------------------------------------------------------
##Prior predictive checks 

#Prior Predictive Check beta0, beta1 ---- FIGURE 3

nsim <- 5000
beta0.prior <- rnorm(nsim, mean=3300, sd=200)
beta1.prior  <- rexp(nsim, rate=1)
slope  <-  beta1.prior
beta0.prior  <- beta0.prior[which((1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:4600
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] - slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure2.1.pdf'), width=6, height=6)
plot(NULL, xlim=c(0,4500), ylim=c(3400,1300), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500), cex.axis=0.9) #X-axis
axis(2, at=seq(3400, 1400, -400))
axis(4, at=BCADtoBP(c(-1400, -1000, -600, -200, 200, 600)), labels=c('1400BC','1000BC','600BC','200BC','200AD','600AD'), cex.axis=0.9)
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=3300, b=-1/0.5, lty=2)
text(x=500, y=1600, label='0.5km/yr')

abline(a=3300, b=-1, lty=2)
text(x=1700, y=1900, label='1km/yr')

abline(a=3300, b=-1/3, lty=2)
text(x=3500, y=2000, label='3km/yr')

abline(a=3300, b=-1/5, lty=2)
text(x=2250, y=2950, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()




