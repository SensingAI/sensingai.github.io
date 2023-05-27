#
# R file to estimate the Lubik-Schorfheide (2007) model
#
rm(list=ls())
library(BMR)
data(BMRLSData)
#
#
#
partomats <- function(parameters){
  psi1     <- parameters[1]
  psi2     <- parameters[2]
  psi3     <- parameters[3]
  #
  rhoR     <- parameters[4]
  alpha    <- parameters[5]
  rSS      <- parameters[6]
  kappa    <- parameters[7]
  tau      <- parameters[8]
  #
  rhoq     <- parameters[9]
  rhoz     <- parameters[10]
  rhoys    <- parameters[11]
  rhopis   <- parameters[12]
  #
  sigmaq   <- parameters[13]^2
  sigmaz   <- parameters[14]^2
  sigmar   <- parameters[15]^2
  sigmays  <- parameters[16]^2
  sigmapis <- parameters[17]^2
  #
  beta     <- exp(-rSS/400)
  tau.alpha   <- tau + alpha*(2 - alpha)*(1 - tau)
  #
  #
  #
  A <- matrix(0,nrow=0,ncol=0)
  B <- matrix(0,nrow=0,ncol=0)
  C <- matrix(0,nrow=0,ncol=0)
  D <- matrix(0,nrow=0,ncol=0)
  #
  F12 <- tau.alpha
  F22 <- beta
  #Order:             y,      pi,         R,         e,        x
  F <- rbind(c(       1,     F12,         0,         0,        0),
             c(       0,     F22,         0,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0))
  #
  G13 <- -tau.alpha
  G21 <- kappa/tau.alpha; G25 <- -kappa/tau.alpha
  G31 <- (1-rhoR)*psi2; G32 <- (1-rhoR)*psi1; G34 <- (1-rhoR)*psi3;
  G <- rbind(c(      -1,       0,       G13,         0,        0),
             c(     G21,      -1,         0,         0,      G25),
             c(     G31,     G32,        -1,       G34,        0),
             c(       0,       1,         0,        -1,        0),
             c(       0,       0,         0,         0,       -1))
  #
  H33 <- rhoR
  H <- rbind(c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,       H33,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0))
  #
  J <- matrix(0,nrow=0,ncol=0)
  K <- matrix(0,nrow=0,ncol=0)
  #
  #Order:             q,      ys,       pis,         z,        r
  L <- rbind(c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        0))
  #
  M11 <- -alpha*tau.alpha*rhoq; M12 <- alpha*(1-alpha)*((1-tau)/tau)*rhoys; M14 <- -rhoz
  M21 <- alpha*(beta*rhoq - 1)
  M41 <- -(1-alpha)
  M52 <- -alpha*(2-alpha)*((1-tau)/tau)
  M <- rbind(c(     M11,     M12,         0,       M14,        0),
             c(     M21,       0,         0,         0,        0),
             c(       0,       0,         0,         0,        1),
             c(     M41,       0,        -1,         0,        0),
             c(       0,     M52,         0,         0,        0))
  #
  N <- rbind(c(    rhoq,       0,         0,         0,        0),
             c(       0,   rhoys,         0,         0,        0),
             c(       0,       0,    rhopis,         0,        0),
             c(       0,       0,         0,      rhoz,        0),
             c(       0,       0,         0,         0,        0))
  #
  shocks <- rbind(c(  sigmaq,       0,         0,         0,        0),
                  c(       0, sigmays,         0,         0,        0),
                  c(       0,       0,  sigmapis,         0,        0),
                  c(       0,       0,         0,    sigmaz,        0),
                  c(       0,       0,         0,         0,   sigmar))
  #
  ObsCons  <- matrix(0,nrow=5,ncol=1)
  MeasErrs <- matrix(0,nrow=5,ncol=5)
  #
  return=list(A=A,B=B,C=C,D=D,F=F,G=G,H=H,J=J,K=K,L=L,M=M,N=N,shocks=shocks,ObsCons=ObsCons,MeasErrs=MeasErrs)
}
#
ObserveMat <- cbind(c(1,0,0,0,0,0,0,0,1,0),
                    c(0,1,0,0,0,0,0,0,0,0),
                    c(0,0,1,0,0,0,0,0,0,0),
                    c(0,0,0,1,0,0,0,0,0,0),
                    c(0,0,0,0,0,1,0,0,0,0))
#
initialvals <- c(1.30,0.23,0.14,
                 0.69,0.11,2.50,0.32,0.31,
                 0.31,0.42,0.97,0.46,
                 1.25,0.84,0.36,1.29,2.00)
#
parnames  <- c( "psi.1","psi.2","psi.3",
                "rhoR","alpha","rSS","kappa","tau",
                "rho.Q","rho.Z","rho.Y","rho.pi",
                "sigma.Q","sigma.Z","sigma.R","sigma.Y","sigma.pi")
#
priorform <- c("Normal","Normal","Normal",
               "Beta","Beta","Normal","Normal","Normal",
               "Beta","Beta","Beta","Beta",
               "Gamma","Gamma","Gamma","Gamma","Gamma")
#
priorpars <- cbind(c(  1.30,   0.23,   0.23,   
                       7,      3,   2.50,  0.32,  0.31,
                       7,      7,     10,     7,
                       2,      2,      2,     2,     3),
                   c(0.20^2, 0.05^2, 0.05^2,   
                     7,      7, 0.15^2, 0.1^2,  0.1^2,
                     7,      7,      4,     7,
                     2,      2,      2,     2,      2))
parbounds <- cbind(c(   NA,    NA,    NA, 
                        0.01,  0.01,    NA,    NA,  NA, 
                        0.01,  0.01,  0.01,  0.01,
                        NA,    NA,    NA,    NA,  NA),
                   c(   NA,    NA,    NA,
                        0.999, 0.999,    NA,    NA,  NA,
                        0.999, 0.999, 0.999, 0.999,
                        NA,    NA,    NA,    NA,  NA))
#
LSest <- EDSGE(LSData,chains=1,cores=1,
               ObserveMat,initialvals,partomats,
               priorform,priorpars,parbounds,parnames,
               optimMethod=c("Nelder-Mead","Nelder-Mead",
                             "Nelder-Mead","Nelder-Mead",
                             "CG"),
               optimLower=NULL,optimUpper=NULL,
               optimControl=list(maxit=10000),
               DSGEIRFs=TRUE,irf.periods=40,
               scalepar=0.25,keep=20000,burnin=30000)
#
parexpre  <- expression(psi[1],psi[2],psi[3],rho[R],alpha,r^{SS},kappa,tau,
                        rho[Q],rho[Z],rho[Y],rho[pi],
                        sigma[Q],sigma[Z],sigma[R],sigma[Y],sigma[pi])
#
modecheck(LSest,gridsize=200,1,FALSE,parnames=parexpre,save=FALSE)
plot(LSest,parnames=parnames,save=FALSE)
#
varnames <- c("Output","Inflation","IntRate","ExchangeRate","PotentialY",
              "ToT","OutputF","InflationF","Technology","IntRateShock")
#
states(LSest,percentiles=c(0.05,0.50,0.95),varnames=varnames)
IRF(LSest,FALSE,varnames=varnames,save=F)
forecast(LSest,5,backdata=5)
#
# DSGE-VAR
#
LSVAR <- DSGEVAR(LSData,chains=1,cores=1,lambda=1,p=4,
                 FALSE,ObserveMat,initialvals,partomats,
                 priorform,priorpars,parbounds,parnames,
                 optimMethod=c("Nelder-Mead","Nelder-Mead",
                               "Nelder-Mead","Nelder-Mead",
                               "CG"),
                 optimLower=NULL,optimUpper=NULL,
                 optimControl=list(maxit=20000),
                 IRFs=TRUE,irf.periods=40,
                 scalepar=0.25,keep=25000,burnin=25000)
#
LSVARInf <- DSGEVAR(LSData,chains=1,cores=1,lambda=Inf,p=4,
                    FALSE,ObserveMat,initialvals,partomats,
                    priorform,priorpars,parbounds,parnames,
                    optimMethod=c("Nelder-Mead","Nelder-Mead",
                                  "Nelder-Mead","CG"),
                    optimLower=NULL,optimUpper=NULL,
                    optimControl=list(maxit=20000),
                    IRFs=TRUE,irf.periods=40,
                    scalepar=0.27,keep=25000,burnin=25000)
#
# To save the graphs, change save=FALSE to save=TRUE in each of the following.
#
modecheck(LSVAR,200,1,FALSE,parnames=parexpre,save=FALSE)
plot(LSVAR,parnames=parexpre,save=FALSE)
#
states(LSVAR,percentiles=c(0.05,0.50,0.95),varnames=varnames)
IRF(LSVAR,varnames=colnames(LSData),comparison=FALSE,save=FALSE)
IRF(LSVAR,varnames=colnames(LSData),comparison=TRUE,save=FALSE)
#
modecheck(LSVARInf,200,1,FALSE,parnames=parexpre,save=FALSE)
plot(LSVARInf,parnames=parexpre,save=FALSE)
#
states(LSVARInf,percentiles=c(0.05,0.50,0.95),varnames=varnames)
IRF(LSVARInf,varnames=colnames(LSData),comparison=FALSE,save=FALSE)
IRF(LSVARInf,varnames=colnames(LSData),comparison=TRUE,save=FALSE)
#
forecast(LSVAR,5,backdata=5)
forecast(LSVARInf,5,backdata=5)
forecast(LSVARInf,5,shocks=FALSE,backdata=5)
#
#
#END