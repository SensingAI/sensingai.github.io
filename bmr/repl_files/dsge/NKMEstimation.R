#
# R file to estimate the basic NK model with BMR
# using Uhlig's method as the solver
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
  A <- matrix(0,nrow=0,ncol=0)
  B <- matrix(0,nrow=0,ncol=0)
  C <- matrix(0,nrow=0,ncol=0)
  D <- matrix(0,nrow=0,ncol=0)
  #
  #Order:            yg        y,        pi,        rn,        i,           n
  F <- rbind(c(      -1,       0,      -eta,         0,        0,           0),
             c(       0,       0,     -beta,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0))
  #
  G <- rbind(c(       1,       0,         0,    -1*eta,      eta,           0),
             c(  -kappa,       0,         1,         0,        0,           0), 
             c(  -phi_y,       0,   -phi_pi,         0,        1,           0),
             c(       0,       0,         0,         1,        0,           0),
             c(       0,       1,         0,         0,        0,  -(1-alpha)),
             c(       1,      -1,         0,         0,        0,           0))
  #
  H <- rbind(c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0),
             c(       0,       0,         0,         0,        0,           0))
  #
  J <- matrix(0,nrow=0,ncol=0)
  K <- matrix(0,nrow=0,ncol=0)
  #
  L <- matrix(c( 0, 0, 
                 0, 0, 
                 0, 0,
                 0, 0, 
                 0, 0,
                 0, 0 ),ncol=2,byrow=T);
  #
  M41 <- -(1/eta)*psi*(rho_a - 1)
  M<- matrix(c(   0,  0,
                  0,  0,
                  0, -1,
                  M41,  0,
                  -1,  0,
                  psi,  0),ncol=2,byrow=T)
  #
  N = matrix(c( rho_a,     0,  
                0, rho_v),nrow=2)
  #
  shocks <- matrix(c(sigmaT,      0,  
                     0, sigmaM),nrow=2)
  #
  ObsCons  <- matrix(0,nrow=2,ncol=1)
  MeasErrs <- matrix(0,nrow=2,ncol=2)
  #
  return=list(A=A,B=B,C=C,D=D,F=F,G=G,H=H,J=J,K=K,L=L,M=M,N=N,shocks=shocks,ObsCons=ObsCons,MeasErrs=MeasErrs)
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
parbounds <- cbind(c(NA,   0.95,  NA, 0.01,  NA,  NA,  NA,  NA, 0.01, 0.01, NA, NA),
                   c(NA,0.99999,  NA,0.999, NA,  NA,  NA,  NA, 0.999,0.999, NA, NA))
#
NKMest <- EDSGE(dsgedata,chains=1,cores=1,
                ObserveMat,initialvals,partomats,
                priorform,priorpars,parbounds,parnames,
                optimMethod=c("Nelder-Mead","CG"),
                optimLower=NULL,optimUpper=NULL,
                optimControl=list(maxit=10000),
                DSGEIRFs=TRUE,irf.periods=40,
                scalepar=0.27,keep=50000,burnin=75000,
                tables=TRUE)
#
parexpre <- expression(alpha,beta,vartheta,theta,eta,phi,phi[pi],phi[Y],
                       rho[A],rho[V],sigma[T],sigma[M])
#
modecheck(NKMest,200,1,FALSE,parnames=parexpre,save=F)
plot(NKMest,parnames=parexpre,save=F)
#
varnames=c("Output Gap","Output","Inflation","Natural Int",
           "Nominal Int","Labour Supply","Technology","MonetaryPolicy")
#
states(NKMest,percentiles=c(0.01,0.50,0.99),varnames=varnames)
IRF(NKMest,FALSE,varnames=varnames,save=F)
forecast(NKMest,backdata=10,save=F)
#
#
#END