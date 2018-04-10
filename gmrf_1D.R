library(TMB)
library(INLA)
setwd("C:/Users/ahav488/Documents/GitHub/gmrf_1D")

x <- seq(1,500,by=1)
n.samp <- 100
sample.x <- sort( sample( 2:(max(x)-1),n.samp ) )

Lambda <- 3
Range <- 30
Kappa <- 2/Range #sqrt(8*nu)/phi, nu=0.5
Tau <- sqrt(gamma(0.5)/(4*pi*Kappa)) #sig2 = 1


mesh.1d <- inla.mesh.1d(sample.x)
spde.1d <- inla.spde2.matern(mesh.1d, alpha = 1)

dyn.load(dynlib("gmrf_1D"))
Data <- list(Y_i = rep(0,n.samp),
             M0 = spde.1d$param.inla$M0,
             M2 = spde.1d$param.inla$M2)
Params <- list(lambda = Lambda, ln_kappa = log(Kappa), 
               ln_tau = log(Tau), omega_i = rep(0,mesh.1d$n))
obj <- MakeADFun(Data, Params, DLL = "gmrf_1D")
sim <- obj$simulate(obj$par)

Data$Y_i <- sim$Y_i
Params <- list(lambda = 0, ln_kappa = 0, ln_tau = 0,
               omega_i = rep(0,mesh.1d$n))
obj <- MakeADFun(Data, Params, random = "omega_i", DLL = "gmrf_1D")
opt <- nlminb(obj$par, obj$fn, obj$gr)

#Check model convergence and fit
opt$message
sdr <- sdreport(obj)
sdr
exp(opt$par[2:3])
report <- obj$report()
report$range
plot(sim$Y_i, exp(report$eta))
plot(sim$Y_i, rpois(n.samp, exp(report$eta)))



#Confirm tmb calculation of Q matches INLA
x <- seq(1,50,by=1)
n.samp <- 10
sample.x <- sort( sample( 2:(max(x)-1),n.samp ) )

Lambda <- 3
Range <- 30
Kappa <- 2/Range #sqrt(8*nu)/phi, nu=0.5
Tau <- sqrt(gamma(0.5)/(4*pi*Kappa)) #sig2 = 1


mesh.1d <- inla.mesh.1d(sample.x)
spde.1d <- inla.spde2.matern(mesh.1d, alpha = 1)
#Equation used in TMB
Q.tmb <- Kappa^2*spde.1d$param.inla$M0 + spde.1d$param.inla$M2 
#Verify with INLA
Q.inla <- inla.spde.precision(spde.1d, theta=c(log(1), log(Kappa)))
Q.tmb;Q.inla
