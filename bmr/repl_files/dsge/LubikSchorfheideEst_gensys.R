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
  G0_13 <- tau.alpha
  G0_16 <- alpha*tau.alpha*rhoq
  G0_17 <- - alpha*(1-alpha)*((1-tau)/tau)*rhoys
  G0_111 <- - tau.alpha
  
  G0_21 <- - kappa/tau.alpha
  G0_25 <- kappa/tau.alpha
  G0_26 <- - alpha*(beta*rhoq - 1)
  
  G0_31 <- - (1-rhoR)*psi2
  G0_32 <- - (1-rhoR)*psi1
  G0_34 <- - (1-rhoR)*psi3
  
  G0_57 <- alpha*(2-alpha)*((1-tau)/tau)
  #                        y,      pi,       R,       e,       x,       q,     y^s,    pi^s,       z,       r,   y_t+1,  pi_t+1,
  Gamma0 <- rbind(c(       1,       0,   G0_13,       0,       0,   G0_16,   G0_17,       0,    rhoz,       0,      -1,  G0_111),
                  c(   G0_21,       1,       0,       0,   G0_25,   G0_26,       0,       0,       0,       0,       0,   -beta),
                  c(   G0_31,   G0_32,       1,   G0_34,       0,       0,       0,       0,       0,      -1,       0,       0),
                  c(       0,      -1,       0,       1,       0, 1-alpha,       0,       1,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       1,       0,   G0_57,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0),
                  c(       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0))
  #
  Gamma1 <- rbind(c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,    rhoR,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,    rhoq,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,   rhoys,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,  rhopis,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,    rhoz,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0),
                  c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1))
  #
  C <- matrix(0,ncol(Gamma0),1)
  #
  Psi <- matrix(0,ncol(Gamma0),5)
  Psi[6:10,] <- diag(5)
  #     
  Pi <- matrix(0,ncol(Gamma0),2)
  Pi[11,1] <- 1
  Pi[12,2] <- 1
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
  return=list(Gamma0=Gamma0,Gamma1=Gamma1,C=C,Psi=Psi,Pi=Pi,shocks=shocks,ObsCons=ObsCons,MeasErrs=MeasErrs)
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
                             "Nelder-Mead","CG"),
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
plot(LSest,parnames=parexpre,save=FALSE)
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
                               "Nelder-Mead","CG"),
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