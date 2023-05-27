#
# An-Schorfheide (2007) Model
#
rm(list=ls())
library(BMR)
#
# Setting parameter values as per An and Schorfheide (2007, Table 2)
#
tau      <- 2.00
kappa    <- 0.15
#
psi1     <- 1.50
psi2     <- 1.00
#
rhoR     <- 0.60
rhog     <- 0.95
rhoz     <- 0.65
#
rA       <- 0.40
piA      <- 4.00
gammaQ   <- 0.50
#
sigmar   <- 0.002
sigmag   <- 0.008
sigmaz   <- 0.0045
#
gamma    <- 1 + (gammaQ/100)
beta     <- 1/( 1 + (rA/400) )
pi       <- 1 + (piA/400)
#
#
#
A <- matrix(0,nrow=0,ncol=0)
B <- matrix(0,nrow=0,ncol=0)
C <- matrix(0,nrow=0,ncol=0)
D <- matrix(0,nrow=0,ncol=0)
#
F12 <- 1/tau
F22 <- beta
#Order:             c,      pi,         y,         R
F <- rbind(c(       1,     F12,         0,         0),
           c(       0,     F22,         0,         0),
           c(       0,       0,         0,         0),
           c(       0,       0,         0,         0))
#
G14 <- -1/tau
G23 <- kappa
G42 <- (1-rhoR)*psi1; G43 <- (1-rhoR)*psi2
G <- rbind(c(      -1,       0,         0,       G14),
           c(       0,      -1,       G23,         0),
           c(       1,       0,        -1,         0),
           c(       0,     G42,       G43,        -1))
#
H44 <- rhoR
H <- rbind(c(       0,       0,         0,         0),
           c(       0,       0,         0,         0),
           c(       0,       0,         0,         0),
           c(       0,       0,         0,       H44))
#
J <- matrix(0,nrow=0,ncol=0)
K <- matrix(0,nrow=0,ncol=0)
#
#Order:             r,       g,         z
L <- rbind(c(       0,       0,         0),
           c(       0,       0,         0),
           c(       0,       0,         0),
           c(       0,       0,         0))
#
M13 <- (1/tau)*rhoz
M22 <- -kappa
M42 <- -(1-rhoR)*psi2
M <- rbind(c(       0,       0,       M13),
           c(       0,     M22,         0),
           c(       0,       1,         0),
           c(       1,     M42,         0))
#
N <- rbind(c(       0,       0,         0),
           c(       0,    rhog,         0),
           c(       0,       0,      rhoz))
#
Sigma <- rbind(c(  sigmar,       0,         0),
               c(       0,  sigmag,         0),
               c(       0,       0,    sigmaz))
#
dsgetest <- uhlig(A,B,C,D,F,G,H,J,K,L,M,N)
#
varnames <- c("Consumption","Inflation","Output","IntRate","IntRateS","GovSpending","Technology")
#
dsgetestirf <- IRF(dsgetest,diag(Sigma),30,varnames=varnames,save=FALSE)
dsgetestsim <- DSGESim(dsgetest,Sigma^2,200,1000,1122,hpfiltered=FALSE)
#
y  <- gammaQ + 100*(dsgetestsim[,3] + dsgetestsim[,7])
pi <- piA + 400*dsgetestsim[,2]
R  <- piA + rA + 4*gammaQ + 400*dsgetestsim[,4]
#
ASData <- cbind(y,pi,R)
#
#save(ASData,file="ASData.RData")
#
#
#END