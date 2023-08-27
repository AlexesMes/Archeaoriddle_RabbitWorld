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
library(raster)
library(terra)
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
win <- win.sf %>% as("Spatial")

#===============================================================================
##Bayesian Quantile Regression
#===============================================================================

## quantreg.R Figures ----

# Bayesian Quantile Regression Model with Measurement Error Plot ---- FIGURE 1
## Compute Fitted Model Confidence Intervals

# rq and median calibrated date
rq.ci <- predict.rq(fit_rq, newdata=data.frame(dist_org=0:620), interval='confidence') #Furthest site from proposed origin, Sq12_1, is 603km

# Bayesian model 
qr.ch1 <- do.call(rbind,quantreg_sample)
post.alpha.quantreg <- qr.ch1[,'alpha']
post.beta.quantreg <- qr.ch1[,'beta']
post.alpha.beta.quantreg  <- data.frame(alpha=post.alpha.quantreg, beta=post.beta.quantreg)
post.theta.quantreg  <- qr.ch1[,grep('theta', colnames(qr.ch1))]
post.theta.med <- apply(post.theta.quantreg, 2, median)
post.quantreg <- apply(post.alpha.beta.quantreg, 1, function(x){x[1]-x[2]*0:620})
post.ci <- t(apply(post.quantreg, 1, quantile, c(0.025,0.5,0.975)))

## Plot and compare
# Transparency color utility function
col.alpha <- function(x,a=1){xx=col2rgb(x)/255;return(rgb(xx[1],xx[2],xx[2],a))}

pdf(file=here('output','figures','figure1.pdf'), width=8.5, height=7)
plot(NULL, xlim=c(0,620), ylim=c(8500,7000), axes=F, xlab='Distance from Sq12_1 Site (in km)', ylab='Cal BP') 
axis(1, at=c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650)) #X-axis
axis(2, at=c(7000, 7250, 7500, 7750, 8000, 8250, 8500)) #Y-axis left
axis(4, at=BCADtoBP(c(-6600, -6400, -6200, -6000, -5800, -5600, -5400, -5200, -5000)), labels=c('6600BC', '6400BC', '6200BC', '6000BC', '5800BC', '5600BC', '5400BC', '5200BC', '5000BC'), cex.axis=0.6) #Y-axis right
points(f_constants$dist_org, f_siteInfo$earliest) #Median Calibrated Date
points(f_constants$dist_org, post.theta.med, pch=20) #Median Posterior theta
for (i in 1:nrow(f_siteInfo))
{
  lines(rep(f_constants$dist_org[i],2), c(f_siteInfo$earliest[i], post.theta.med[i]),lty=2)
}
lines(0:620, rq.ci[,1], lty=1, lwd=2, col='blue') #Quantile Regression on Median Dates
polygon(x=c(0:620, 620:0), c(rq.ci[,2], rev(rq.ci[,3])), col=col.alpha('lightblue', 0.4), border=NA)

lines(0:620, post.ci[,2], lty=1, lwd=2, col='indianred') #Bayesian Quantile Regression with Measurement Error
polygon(x=c(0:620, 620:0), c(post.ci[,1], rev(post.ci[,3])), col=col.alpha('indianred', 0.2), border=NA)


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
beta0.prior <- rnorm(nsim, mean=8000, sd=200)
beta1.prior  <- rexp(nsim, rate=1)
slope  <-  beta1.prior
beta0.prior  <- beta0.prior[which((1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:620
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] - slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure3.pdf'), width=8, height=8)
plot(NULL, xlim=c(0,620), ylim=c(8500,7000), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650)) #X-axis
axis(2, at=c(7000, 7250, 7500, 7750, 8000, 8250, 8500)) #Y-axis left
axis(4, at=BCADtoBP(c(-6600, -6400, -6200, -6000, -5800, -5600, -5400, -5200, -5000)), labels=c('6600BC', '6400BC', '6200BC', '6000BC', '5800BC', '5600BC', '5400BC', '5200BC', '5000BC'), cex.axis=0.6) #Y-axis right
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=8000, b=-1/0.5, lty=2)
text(x=200, y=7500, label='0.5km/yr')

abline(a=8000, b=-1, lty=2)
text(x=270, y=7670, label='1km/yr')

abline(a=8000, b=-1/3, lty=2)
text(x=300, y=7850, label='3km/yr')

abline(a=8000, b=-1/5, lty=2)
text(x=350, y=7990, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()



#===============================================================================
###Gaussian Process Quantile Regression
#===============================================================================

##Effect of parameter variation in covariance matrix
# Variance Covariance Function Parameters ---- FIGURE 4
pdf(here('output','figures','figure4.pdf'), width=10, height=5)

layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))

plot(NULL, xlim=c(0,500), ylim=c(0.00, 0.05), xlab=TeX('$d_{i,j}$'), ylab=TeX('$cov(i,j)$'))

#Default parameter values
etasq_def = 0.05
rho_def = 100
#Varying parameter values
rho = c(25, 50, 75, 100, 125)
etasq = c(0.025, 0.05, 0.075, 0.1)

cl <- rainbow(5) #Assigning random colours
for (i in 1:length(rho))
{
  curve(etasq_def*exp(-0.5*(x/rho[i])^2), from=0, to=500, col=cl[i], add = TRUE)
}
abline(h=0.05,lty=2)
lines(x=c(rho_def, rho_def), y=c(-0.01, etasq_def*exp(-0.5)),lty=2)
lines(x=c(-20, rho_def), y=c(etasq_def*exp(-0.5), etasq_def*exp(-0.5)),lty=2)

text(400, 0.046, label=TeX('$\\eta^2 = 0.05$'))
text(75, etasq_def*exp(-0.5)-0.004, label=TeX('$\\eta^2e^{-0.5}$'))
text(27, 0.006, label=TeX('$\\rho = 25$'), cex=0.7)
text(82, 0.006, label=TeX('$\\rho = 50$'), cex=0.7)
text(150, 0.006, label=TeX('$\\rho = 75$'), cex=0.7)
text(225, 0.006, label=TeX('$\\rho = 100$'), cex=0.7)
text(305, 0.006, label=TeX('$\\rho = 125$'), cex=0.7)

plot(NULL, xlim=c(0,500), ylim=c(0.00, 0.1), xlab=TeX('$d_{i,j}$'), ylab=TeX('$cov(i,j)$'), cex.sub = 0.6) #sub=TeX("The quadratic exponential model and its relationship to the parameters $\\rho$, $\\eta$, and $d_{i,j}$, assuming ${i \\neq j}$.")
for (i in 1:length(etasq))
{
  curve(etasq[i]*exp(-0.5*(x/rho_def)^2), from=0, to=500, col=cl[i], add = TRUE)
}
text(450, 0.096, label=TeX('$\\rho = 100$'))
text(40, 0.015, label=TeX('$\\eta^2 = 0.025$'), cex=0.7)
text(40, 0.04, label=TeX('$\\eta^2 = 0.05$'), cex=0.7)
text(40, 0.06, label=TeX('$\\eta^2 = 0.075$'), cex=0.7)
text(40, 0.09, label=TeX('$\\eta^2 = 0.1$'), cex=0.7)

dev.off()


#-------------------------------------------------------------------------------
#Impact of rho and etasq on variability in dispersal rate ---- FIGURE 5
source(here('src','gpqrSim.R'))

# Fixed parameters
n = 300
origin.point = c(11.40, 5.48)
beta0 = 8000
beta1 = 0.5
sigma = 100
seed = 144231

# Sweep parameters
etasq = c(0.01, 0.05, 0.15)
rho = c(20, 75, 150)
cov_param = expand.grid(etasq=etasq, rho=rho)
tmp  <- vector('list',length=nrow(cov_param))
out  <- vector('list',length=nrow(cov_param))


# Simulate sites and plot
for (i in 1:nrow(cov_param))
{
  tmp[[i]]  <- gpqrSim(win = win,
                       n = n,
                       beta0 = beta0,
                       beta1 = beta1,
                       sigma = sigma,
                       origin.point = origin.point,
                       etasq = cov_param$etasq[i],
                       rho =cov_param$rho[i],
                       seed = seed)
  
  out[[i]] <- ggplot() +
    geom_sf(data=win.sf,aes(), fill='grey66', show.legend=FALSE, lwd=0) +
    geom_sf(data=tmp[[i]], mapping = aes(fill=rate), pch=21, col='darkgrey', size=1.5) +
    xlim(-4.5,1.5) + #Center the frame on RabbitWorld
    ylim(-1.5,4.5) +
    labs(title=TeX(sprintf("$\\eta^2 = %g \\, \\rho = %g$", cov_param$etasq[i], cov_param$rho[i])), fill='Dispersal Rate \n (km/yr)') +
    scale_fill_viridis(option = 'turbo', limits = c(0.7, 4), oob = scales::squish) +
    theme(plot.title = element_text(hjust = 0.5, size=11), panel.background = element_rect(fill='lightblue'), panel.grid.major = element_line(size = 0.1), legend.position=c(0.2, 0.2), legend.text = element_text(size=7), legend.key.width= unit(0.1, 'in'), legend.key.size = unit(0.08, "in"), legend.background=element_rect(fill = alpha("white", 0.5)), legend.title=element_text(size=7), axis.text=element_blank(), axis.ticks=element_blank(), plot.margin = unit(c(0,0,0,0), "in"))
  
}

pdf(file=here('output','figures','figure5.pdf'), width=8, height=8)
grid.arrange(out[[1]], out[[2]], out[[3]], out[[1]], out[[4]], out[[7]], nrow=2, ncol=3)
dev.off()

#-------------------------------------------------------------------------------
# Relationship between slope and rate of dispersal ---- FIGURE 6
pdf(here('output','figures','figure6.pdf'), width=5, height=5)
par(mar=c(5,6,2,2))
curve(-1/x, from=-2.5, to=-0.05, xlab=TeX('$Slope \\, s - \\beta_1$'), ylab=TeX('$Speed \\, \\frac{-1}{s-\\beta_1}$'), axes=TRUE, sub=TeX("The relationship between regression slope and its negative reciprocal, i.e. the rate of dispersal."), cex.sub = 0.6)
dev.off()


#-------------------------------------------------------------------------------
##Prior predictive checks 

#Prior Predictive Check beta0, beta1, s ---- FIGURE 7
nsim <- 5000
beta0.prior <- rnorm(nsim, mean=8000, sd=200)
beta1.prior  <- rexp(nsim, rate=2)
s.prior  <- rnorm(nsim, mean=0, sd=sqrt(rexp(nsim, rate=20)))
slope  <- s.prior - beta1.prior
beta0.prior  <- beta0.prior[which((-1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((-1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:620
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] + slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure7.pdf'), width=8, height=8)
plot(NULL, xlim=c(0,620), ylim=c(8500,7000), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650)) #X-axis
axis(2, at=c(7000, 7250, 7500, 7750, 8000, 8250, 8500)) #Y-axis left
axis(4, at=BCADtoBP(c(-6600, -6400, -6200, -6000, -5800, -5600, -5400, -5200, -5000)), labels=c('6600BC', '6400BC', '6200BC', '6000BC', '5800BC', '5600BC', '5400BC', '5200BC', '5000BC'), cex.axis=0.6) #Y-axis right
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=8000, b=-1/0.5, lty=2)
text(x=200, y=7500, label='0.5km/yr')

abline(a=8000, b=-1, lty=2)
text(x=270, y=7670, label='1km/yr')

abline(a=8000, b=-1/3, lty=2)
text(x=300, y=7850, label='3km/yr')

abline(a=8000, b=-1/5, lty=2)
text(x=350, y=7990, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

#-------------------------------------------------------------------------------
# Prior Predictive Check etasq and rho ---- FIGURE 8
nsim  <- 1000
etasq.prior  <- rexp(nsim, 20)
rho.prior  <- rtgamma(nsim,10,(10-1)/100, 1, 620)
cov.mat = matrix(NA, nrow=nsim, ncol=length(0:700))
for (i in 1:nsim)
{
  cov.mat[i,] = etasq.prior[i]*exp(-0.5*(0:700/rho.prior[i])^2)
}

pdf(file=here('output', 'figures', 'figure8.pdf'), width=6, height=5)
plot(NULL, xlab='Distance (km)', ylab='Covariance', xlim=c(0,700), ylim=c(0,0.2))
polygon(c(0:700, 700:0), c(apply(cov.mat, 2, quantile, 0.025), rev(apply(cov.mat, 2, quantile,0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(c(0:700, 700:0), c(apply(cov.mat, 2, quantile, 0.5), rev(apply(cov.mat, 2, quantile,0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

legend('topright', legend=c('50% percentile range','95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

#===============================================================================
## qpqr_tau90.R and qpqr_tau99.R Figures ----

# Load data
load(here('output','gpqr_tau90.RData'))
load(here('output','gpqr_tau99.RData'))

#-------------------------------------------------------------------------------
# Posterior Mean of dispersal rate deviations ---- FIGURE 8

post.gpqr.tau90  <- do.call(rbind,gpqr_tau90)
post.gpqr.tau99  <- do.call(rbind,gpqr_tau99)

nmcmc  <- nrow(post.gpqr.tau90) #number of MCMC samples (same for tau90 and tau99)
post.s.tau90  <- post.gpqr.tau90[,grep('s\\[',colnames(post.gpqr.tau90))]
post.s.tau99  <- post.gpqr.tau99[,grep('s\\[',colnames(post.gpqr.tau99))]

#post.arrival <- matrix(NA,nmcmc,constants$N.sites)
post.rate.tau90 <- post.rate.tau99  <-  matrix(NA,nmcmc,f_constants$n_sites)
for (i in 1:nmcmc)
{
  post.rate.tau90[i,] = -1 / (post.s.tau90[i,]-post.gpqr.tau90[i,'beta1'])
  post.rate.tau99[i,] = -1 / (post.s.tau99[i,]-post.gpqr.tau99[i,'beta1'])
}

f_sites@data$s.m.tau90 <- apply(post.s.tau90,2,median)
f_sites@data$s.m.tau99 <- apply(post.s.tau99,2,median)
f_sites@data$rate.m.tau90  <- apply(post.rate.tau90,2,median)
f_sites@data$rate.m.tau99  <- apply(post.rate.tau99,2,median)

sites.sf <- as(f_sites, 'sf')

f2a <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='A'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win.sf, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=s.m.tau90), pch=21, col='black', size=2) +
  xlim(-4.5,1.5) + #Center the frame on RabbitWorld
  ylim(-1.5,4.5) +
  labs(title=TeX(r"(Posterior Median of s with $\tau = 0.90$)"), fill='s') +
  scale_fill_gradient2(low='blue', high='red', mid='white') +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10),
        panel.background = element_rect(fill='lightblue'),
        panel.grid.major = element_line(size = 0.1),
        legend.position=c(0.13,0.12),
        legend.text = element_text(size=6),
        legend.key.width= unit(0.1, 'in'),
        legend.key.size = unit(0.08, "in"),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.title=element_text(size=7),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

f2b <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='B'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win.sf, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=s.m.tau99), pch=21, col='black', size=2) +
  xlim(-4.5,1.5) + #Center the frame on RabbitWorld
  ylim(-1.5,4.5) +
  labs(title=TeX(r"(Posterior Median of s with $\tau = 0.99$)"), fill='s') +
  scale_fill_gradient2(low='blue', high='red', mid='white') +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10),
        panel.background = element_rect(fill='lightblue'),
        panel.grid.major = element_line(size = 0.1),
        legend.position=c(0.13,0.12),
        legend.text = element_text(size=6),
        legend.key.width= unit(0.1, 'in'),
        legend.key.size = unit(0.08, "in"),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.title=element_text(size=7),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

f2c <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='C'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win.sf, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=rate.m.tau90), pch=21, col='black', size=2) +
  xlim(-4.5,1.5) + #Center the frame on RabbitWorld
  ylim(-1.5,4.5) +
  labs(title=TeX(r"(Posterior median of dispersal rate with $\tau = 0.90$)"), fill='Dispersal Rate (km/year)') +
  scale_fill_viridis(option="turbo", limits=c(0,4)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10),
        panel.background = element_rect(fill='lightblue'),
        panel.grid.major = element_line(size = 0.1),
        legend.position=c(0.18,0.12),
        legend.text = element_text(size=6),
        legend.key.width= unit(0.1, 'in'),
        legend.key.size = unit(0.08, "in"),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.title=element_text(size=6),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

f2d <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='D'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win.sf, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=rate.m.tau99), pch=21, col='black', size=2) +
  xlim(-4.5,1.5) + #Center the frame on RabbitWorld
  ylim(-1.5,4.5) +
  labs(title=TeX(r"(Posterior median of dispersal rate with $\tau = 0.99$)"), fill='Dispersal Rate (km/year)') +
  scale_fill_viridis(option="turbo", limits=c(0,4)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10),
        panel.background = element_rect(fill='lightblue'),
        panel.grid.major = element_line(size = 0.1),
        legend.position=c(0.18,0.12),
        legend.text = element_text(size=6),
        legend.key.width= unit(0.1, 'in'),
        legend.key.size = unit(0.08, "in"),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.title=element_text(size=6),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


pdf(file=here('output','figures','figure8.pdf'), width=7, height=7)
grid.arrange(f2a, f2b, f2c, f2d, ncol=2, padding=0)
dev.off()




#-------------------------------------------------------------------------------
# Traceplot of beta0, beta1, rhosq, and etasq for tau=0.9) ---- FIGURE 9

pdf(file=here('output','figures','figure9.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau90[,'beta0'], main=TeX('$\\beta_0$'), smooth=TRUE)
traceplot(gpqr_tau90[,'beta1'], main=TeX('$\\beta_1$'), smooth=TRUE)
traceplot(gpqr_tau90[,'rho'], main=TeX('$\\rho$'), smooth=TRUE)
traceplot(gpqr_tau90[,'etasq'], main=TeX('$\\eta^2$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Traceplot of beta0, beta1, rhosq, and etasq for tau=0.99 ---- FIGURE 10

pdf(file=here('output','figures','figure10.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau99[,'beta0'], main=TeX('$\\beta_0$'), smooth=TRUE)
traceplot(gpqr_tau99[,'beta1'], main=TeX('$\\beta_1$'), smooth=TRUE)
traceplot(gpqr_tau99[,'rho'], main=TeX('$\\rho$'), smooth=TRUE)
traceplot(gpqr_tau99[,'etasq'], main=TeX('$\\eta^2$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.9 ---- FIGURE 11

gpqr.tau90.comb  <- do.call(rbind, gpqr_tau90)

pdf(file=here('output','figures','figure11.pdf'), width=8, height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau90.comb[,'beta0'], main=TeX('$\\beta_0$'), xlab='Cal BP', ylab='')
postHPDplot(gpqr.tau90.comb[,'beta1'], main=TeX('$\\beta_1$'), xlab='', ylab='')
postHPDplot(gpqr.tau90.comb[,'rho'], main=TeX('$\\rho$'), xlab='km', ylab='')
postHPDplot(gpqr.tau90.comb[,'etasq'], main=TeX('$\\eta^2$'), xlab='', ylab='')
dev.off()

#-------------------------------------------------------------------------------
# Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.99 ---- FIGURE 12

gpqr.tau99.comb  <- do.call(rbind, gpqr_tau99)

pdf(file=here('output','figures','figure12.pdf'), width=8, height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau99.comb[,'beta0'], main=TeX('$\\beta_0$'), xlab='Cal BP', ylab='')
postHPDplot(gpqr.tau99.comb[,'beta1'], main=TeX('$\\beta_1$'), xlab='', ylab='')
postHPDplot(gpqr.tau99.comb[,'rho'], main=TeX('$\\rho$'), xlab='km', ylab='')
postHPDplot(gpqr.tau99.comb[,'etasq'], main=TeX('$\\eta^2$'), xlab='', ylab='')
dev.off()


#===============================================================================
###Bayesian Hierarchical Phase Model
#===============================================================================
##Bayesian Hierarchical Phase Models without constraints

#Load Data ----
load(here("output", "phase_model_a.RData"))

#-------------------------------------------------------------------------------
# Prior Predictive check for duration parameter, delta ---- FIGURE 13
nsim  <- 5000

set.seed(123)

gamma1  <- runif(nsim,1,20)
gamma2  <- rtruncnorm(nsim, mean=200, sd=100, 1, 500)
delta.mat = matrix(NA, ncol=1000, nrow=nsim) #Initialise

for (i in 1:nsim) {
  delta.mat[i,] = dgamma(1:1000, gamma1[i], (gamma1[i]-1)/gamma2[i])
}

pdf(file=here('output','figures','figure13.pdf'), height=6, width=6)

plot(NULL,xlab=TeX('$\\delta$'),ylab='Probability Density',xlim=c(1,1000),ylim=c(0,0.02))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.025), rev(apply(delta.mat,2,quantile,prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.25), rev(apply(delta.mat,2,quantile,prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))
legend('topright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))

dev.off()


#-------------------------------------------------------------------------------
# Marginal Posterior Distribution of nu, model a ---- FIGURE 14

out.comb.unif.model.a  <- do.call(rbind, out_unif_model_a)
post.nu.model.a  <- out.comb.unif.model.a[,paste0('a[',1:25,']')] %>%  round() #100 sq areas
model.a.long  <- data.frame(value = as.numeric(post.nu.model.a),
                            area = factor(rep(1:25, each=nrow(post.nu.model.a)), levels=paste0(1:25), ordered=TRUE))

#Plot
pdf(file=here('output','figures','figure14.pdf'), height=10, width=8)

ggplot(model.a.long, aes(x = value, y = area, fill='lightblue')) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(9000,7000), breaks=BCADtoBP(c(-7000, -6800, -6600, -6400, -6200, -6000, -5800, -5600, -5400, -5200)), labels=c('7000BC', '6800BC', '6600BC', '6400BC', '6200BC', '6000BC', '5800BC', '5600BC', '5400BC', '5200BC')) +
  scale_fill_manual(values='lightblue') +
  theme(legend.position = "none") +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$')))

dev.off()


#-------------------------------------------------------------------------------
# Probability Matrix of nu, model a ---- FIGURE 15
source(here('src','orderPPlot.R'))
post.nu.model.a_rel  <- out.comb.unif.model.a[,paste0('a[',1:25,']')] %>%  round() #Keep relevant sq areas

pdf(file=here('output','figures','figure15.pdf'), width=10, height=10.5)
orderPPlot(post.nu.model.a_rel, name.vec=paste("Area", 1:25))
dev.off()


#-------------------------------------------------------------------------------
## Plot sq areas with median arrival times ---- FIGURE 16

#Extract arrival times for model A
out.comb.unif.modela  <- do.call(rbind, out_unif_model_a)
post.nu.modela  <- out.comb.unif.modela[,paste0('a[',1:25,']')]  %>% round() 
hpdi.modela  <- apply(post.nu.modela, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.modela  <- apply(post.nu.modela, 2, median)
hi90_modA  <- hpdi.modela[1,]
lo90_modA  <- hpdi.modela[2,]

area_sites = c(1,2,4,10,13,14,16,18,24) #areas which contain sites
ocean_boundary = st_difference(sq_grid$geometry, as(gUnaryUnion(sample_win_sp), "sf")$geometry)

median_sq_dates_modA <- sq_grid %>% 
  mutate(median_date = as.factor(med.modela),
         hpdi_high = as.factor(hi90_modA),
         hpdi_low = as.factor(lo90_modA),
         contains_sites = as.factor(case_when(area_id %in% area_sites ~ 1, area_id %!in% area_sites ~ 0)),
         median_date_sq_with_sites = case_when(area_id %in% area_sites ~ median_date, area_id %!in% area_sites ~ NA)) 


#Plot
#-----MODEL A
modA <- ggplot(data = median_sq_dates_modA) +
  geom_sf(aes(fill = median_date, alpha=contains_sites)) + #sq grid
  geom_sf(data = ocean_boundary, color = "lightblue", fill="lightblue") + #block out ocean
  geom_sf(data=sq_grid, fill=NA)+ #outline grid squares
  scale_fill_manual(values = wes_palette(43, name = "GrandBudapest1", type = "continuous"), name = "")+
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = ifelse(is.na(median_date_sq_with_sites), "", paste0(median_date_sq_with_sites, "BP"))), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
  scale_alpha_manual(values=c(0.5, 1)) +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


modAHPDIlow <- ggplot(data = median_sq_dates_modA) +
  geom_sf(aes(fill = hpdi_low, alpha=contains_sites)) +
  geom_sf(data = ocean_boundary, color = "lightblue", fill="lightblue") + #block out ocean
  geom_sf(data=sq_grid, fill=NA)+ #outline grid squares
  ggtitle('90% HPDI (low)') +
  scale_fill_manual(values = wes_palette(43, name = "GrandBudapest1", type = "continuous"), name = "") +
  scale_alpha_manual(values=c(0.5, 1)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")

modAHPDIhigh <- ggplot(data = median_sq_dates_modA) +
  geom_sf(aes(fill = hpdi_high, alpha=contains_sites)) +
  geom_sf(data = ocean_boundary, color = "lightblue", fill="lightblue") + #block out ocean
  geom_sf(data=sq_grid, fill=NA)+ #outline grid squares
  ggtitle('90% HPDI (high)') +
  scale_fill_manual(values = wes_palette(43, name = "GrandBudapest1", type = "continuous"), name = "") +
  scale_alpha_manual(values=c(0.5, 1)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")


A <- cowplot::ggdraw() +
  draw_plot(modA) +
  draw_plot(modAHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modAHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#Output
pdf(file=here('output','figures','map_figure_arrival.pdf'), width=15, height=8)
grid.arrange(A, ncol=1, padding=0)
dev.off()


#-------------------------------------------------------------------------------
## Plot cumulative friction ---- FIGURE 17

#Load Data ----
load(here("output", "cumulative_friction.RData"))


pdf(file=here('output','figures','map_cumulative_friction_100.pdf'), width=15, height=8)

##Plot square areas coloured by site frequency, fitness values, elevation or cumulative friction
ggplot(data = sq_grid_conc) +
  geom_sf(alpha=1, aes(fill=cf)) + #Alternatively: fill=sq_fitness_value #fill=sq_height_value
  labs(x = "Longitude", y = "Latitude", fill='Cumulative Friction')  + 
  geom_sf(data = ocean_boundary, color = "lightblue", fill="lightblue") + #block out ocean
  geom_sf(data = as(f_sites, 'sf'), size=3, alpha=0.5, aes(colour="purple")) + #sites
  guides(colour = "none") + #don't plot site legend
  geom_sf(data=sq_grid, fill=NA)+ #outline grid squares
  geom_sf_label(aes(label = area_id), label.size  = NA, alpha = 0.4, size=3.5) + #sq grid labels
  scale_fill_viridis() +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"))

dev.off()





