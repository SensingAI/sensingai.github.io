#
# R file to estimate a DSGE-VAR for the basic NK model
# using gensys as the solver
#
rm(list=ls())
library(BMR)
# Load and demean the data
data(BMRVARData)
dsgedata <- USMacroData[24:211,-c(1,3)]
dsgedata <- as.matrix(dsgedata)
for(i in 1:2){
  dsgedata[,i] <- dsgedata[,i] - mean(dsgedata[,i])
}
#
#
#
partomats <- function(parameters){
  alpha    <- parameters[1]
  beta     <- parameters[2]
  vartheta <- parameters[3]
  theta    <- parameters[4]
  #
  eta      <- parameters[5]
  phi      <- parameters[6]
  phi_pi   <- parameters[7]
  phi_y    <- parameters[8]
  rho_a    <- parameters[9]
  rho_v    <- parameters[10]
  #
  sigmaT   <- (parameters[11])^2
  sigmaM   <- (parameters[12])^2
  #
  BigTheta <- (1-alpha)/(1-alpha+alpha*vartheta)
  kappa <- (((1-theta)*(1-beta*theta))/(theta))*BigTheta*((1/eta)+((phi+alpha)/(1-alpha)))
  #
  psi <- (eta*(1+phi))/(1-alpha+eta*(phi + alpha))
  #
  G0_47 <- (1/eta)*psi*(rho_a - 1)
  #Order:                 yg,       y,      pi,      rn,       i,       n,       a,       v,  yg_t+1,  pi_t+1
  Gamma0 <- rbind(c(      -1,       0,       0,     eta,  -eta/4,       0,       0,       0,       1,   eta/4),
                  c(   kappa,       0,    -1/4,       0,       0,       0,       0,       0,       0,  beta/4),
                  c(   phi_y,       0,phi_pi/4,       0,    -1/4,       0,       0,       1,       0,       0),
                  c(       0,       0,       0,      -1,       0,       0,   G0_47,       0,       0,       0),
                  c(       0,      -1,       0,       0,       0, 1-alpha,       1,       0,       0,       0),
                  c(      -1,       1,       0,       0,       0,       0,    -psi,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       1,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       1,       0,       0),
                  c(       1,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       1,       0,       0,       0,       0,       0,       0,       0))
  #
  Gamma1 <- rbind(c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,   rho_a,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,   rho_v,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       1,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       1))
  #
  C <- matrix(0,10,1)
  #
  Psi <- matrix(0,10,2)
  Psi[7,1] <- 1
  Psi[8,2] <- 1
  #     
  Pi <- matrix(0,10,2)
  Pi[9,1] <- 1
  Pi[10,2] <- 1
  #
  shocks <- matrix(c(sigmaT,      0,  
                     0, sigmaM),nrow=2)
  #
  ObsCons  <- matrix(0,nrow=2,ncol=1)
  MeasErrs <- matrix(0,nrow=2,ncol=2)
  #
  return=list(Gamma0=Gamma0,Gamma1=Gamma1,C=C,Psi=Psi,Pi=Pi,shocks=shocks,ObsCons=ObsCons,MeasErrs=MeasErrs)
}
#
ObserveMat <- cbind(c(0,0,4,0,0,0,0,0),
                    c(0,0,0,0,4,0,0,0))
#
initialvals <- c(.33,0.97,6,0.6667,1,1,1.5,0.5/4,0.90,0.5,1,0.25)
parnames  <- c( "Alpha","Beta","Vartheta","Theta","Eta","Phi","Phi.Pi","Phi.Y","Rho.A","Rho.V","Sigma.T","Sigma.M")
priorform <- c("Normal","Beta","Normal","Beta","Normal","Normal","Normal","Normal","Beta","Beta","Gamma","IGamma")
priorpars <- cbind(c(  0.33, 20,   5,   3,     1,     1,   1.5,   0.7,   3,   4,  2,  2),
                   c(0.05^2,  2,   1,   2, 0.1^2, 0.3^2, 0.1^2, 0.1^2,   2, 1.5,  2,  1))
parbounds <- cbind(c(NA, 0.95,  NA,0.01,  NA,  NA,  NA,  NA, 0.01, 0.01, NA, NA),
                   c(NA,0.99999,  NA,0.999, NA,  NA,  NA,  NA,0.999,0.999, NA, NA))
#
NKMVAR <- DSGEVAR(dsgedata,chains=1,cores=1,lambda=1,p=4,
                  FALSE,ObserveMat,initialvals,partomats,
                  priorform,priorpars,parbounds,parnames,
                  optimMethod=c("Nelder-Mead","CG"),
                  optimLower=NULL,optimUpper=NULL,
                  optimControl=list(maxit=20000,reltol=(10^(-12))),
                  IRFs=TRUE,irf.periods=5,
                  scalepar=0.28,keep=25000,burnin=25000)
#
NKMVARInf <- DSGEVAR(dsgedata,chains=1,cores=1,lambda=Inf,p=4,
                     FALSE,ObserveMat,initialvals,partomats,
                     priorform,priorpars,parbounds,parnames,
                     optimMethod=c("Nelder-Mead","CG"),
                     optimLower=NULL,optimUpper=NULL,
                     optimControl=list(maxit=20000,reltol=(10^(-12))),
                     IRFs=TRUE,irf.periods=5,
                     scalepar=0.3,keep=25000,burnin=25000)
#
# To save the graphs, change save=FALSE to save=TRUE in each of the following.
#
parexpre <- expression(alpha,beta,vartheta,theta,eta,phi,phi[pi],phi[Y],
                       rho[A],rho[V],sigma[T],sigma[M])
#
modecheck(NKMVAR,200,1,FALSE,parnames=parexpre,save=FALSE)
plot(NKMVAR,parnames=parexpre,save=FALSE)
#
IRF(NKMVAR,varnames=names(USMacroData)[c(2,4)],save=FALSE)
IRF(NKMVAR,varnames=names(USMacroData)[c(2,4)],comparison=FALSE,save=FALSE)
#
modecheck(NKMVARInf,200,1,FALSE,parnames=parexpre,save=FALSE)
plot(NKMVARInf,parnames=parexpre,save=FALSE)
#
IRF(NKMVARInf,varnames=names(USMacroData)[c(2,4)],save=FALSE)
IRF(NKMVARInf,varnames=names(USMacroData)[c(2,4)],comparison=FALSE,save=FALSE)
#
forecast(NKMVAR,periods=5,backdata=5,save=FALSE)
forecast(NKMVAR,periods=5,shocks=FALSE,backdata=5,save=FALSE)
forecast(NKMVARInf,periods=5,backdata=5,save=FALSE)
#
#
# END