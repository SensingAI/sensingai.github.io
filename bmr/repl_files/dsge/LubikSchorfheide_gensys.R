#
# Lubik-Schorfheide (2007) Model
#
rm(list=ls())
library(BMR)
#
# Setting parameter values as per Lubik and Schorfheide (2007, Table 3, posterior means)
#
psi1     <- 1.30
psi2     <- 0.23
psi3     <- 0.14
#
rhoR     <- 0.69
alpha    <- 0.11
rSS      <- 2.52
kappa    <- 0.32
tau      <- 0.31
#
rhoq     <- 0.31
rhoz     <- 0.42
rhoys    <- 0.97
rhopis   <- 0.46
#
sigmar   <- 0.36
sigmaq   <- 1.25
sigmaz   <- 0.84
sigmays  <- 1.29
sigmapis <- 2.00
#
beta      <- exp(-rSS/400)
tau.alpha <- tau + alpha*(2 - alpha)*(1 - tau)
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
Sigma <- rbind(c(  sigmaq,       0,         0,         0,        0),
               c(       0, sigmays,         0,         0,        0),
               c(       0,       0,  sigmapis,         0,        0),
               c(       0,       0,         0,    sigmaz,        0),
               c(       0,       0,         0,         0,   sigmar))
#
dsgetest <- gensys(Gamma0,Gamma1,C,Psi,Pi)
#
varnames <- c("Output","Inflation","IntRate","ExchangeRate","PotentialY",
              "ToT","OutputF","InflationF","Technology","IntRateShock")
#
dsgetestirf <- IRF(dsgetest,diag(Sigma),30,varnames=varnames,save=FALSE)
dsgetestsim <- DSGESim(dsgetest,Sigma^2,200,1000,1122,hpfiltered=FALSE)
#
y  <- dsgetestsim[,1] + dsgetestsim[,9]
pi <- dsgetestsim[,2]
R  <- dsgetestsim[,3]
e  <- dsgetestsim[,4]
q  <- dsgetestsim[,6]
#
LSData <- cbind(y,pi,R,e,q)
#
#save(LSData,file="LSData.RData")
#
#
#END