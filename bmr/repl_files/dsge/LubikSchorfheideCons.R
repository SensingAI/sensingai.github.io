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
# Values taken from Anders Warne's YADA LS Example:
piA      <- 1.95
gammaQ   <- 0.55
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
Sigma <- rbind(c(  sigmaq,       0,         0,         0,        0),
               c(       0, sigmays,         0,         0,        0),
               c(       0,       0,  sigmapis,         0,        0),
               c(       0,       0,         0,    sigmaz,        0),
               c(       0,       0,         0,         0,   sigmar))
#
dsgetest <- uhlig(A,B,C,D,F,G,H,J,K,L,M,N)
#
varnames <- c("Output","Inflation","IntRate","ExchangeRate","PotentialY",
              "ToT","OutputF","InflationF","Technology","IntRateShock")
#
dsgetestirf <- IRF(dsgetest,diag(Sigma),30,varnames=varnames,save=FALSE)
dsgetestsim <- DSGESim(dsgetest,Sigma^2,200,1000,1122,hpfiltered=FALSE)
#
y  <- gammaQ + dsgetestsim[,1] + dsgetestsim[,9]
pi <- piA + dsgetestsim[,2]
R  <- piA + rSS + 4*gammaQ + dsgetestsim[,3]
e  <- dsgetestsim[,4]
q  <- dsgetestsim[,6]
#
LSDataCons <- cbind(y,pi,R,e,q)
#
#save(LSDataCons,file="LSDataCons.RData")
#
#
#END