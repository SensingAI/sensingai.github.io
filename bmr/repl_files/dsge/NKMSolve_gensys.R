#
# Solve the NK model with gensys
#
rm(list=ls())
library(BMR)
#
alpha    <- 0.33
vartheta <- 6
beta     <- 0.99
theta    <- 0.6667
#
eta    <- 1               
phi    <- 1                  
phi_pi <- 1.5             
phi_y  <- 0.5/4
rho_a  <- 0.90
rho_v  <- 0.5
#
BigTheta <- (1-alpha)/(1-alpha+alpha*vartheta)
kappa    <- (((1-theta)*(1-beta*theta))/(theta))*BigTheta*((1/eta)+((phi+alpha)/(1-alpha)))
psi      <- (eta*(1+phi))/(1-alpha+eta*(phi + alpha))
#
sigmaTech     <- 1
sigmaMonetary <- 0.25
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
Sigma <- rbind(c(sigmaTech,             0),  
               c(        0, sigmaMonetary))
#
dsgetest <- gensys(Gamma0,Gamma1,C,Psi,Pi)
#
varnames <- c("Output Gap","Output","Inflation","Natural Int",
              "Nominal Int","Labour Supply","Technology","MonetaryPolicy")
#
dsgetestirf <- IRF(dsgetest,diag(Sigma),12,varnames=varnames,save=FALSE)
dsgetestsim <- DSGESim(dsgetest,Sigma^2,800,200,1122,FALSE)
#
#
#END