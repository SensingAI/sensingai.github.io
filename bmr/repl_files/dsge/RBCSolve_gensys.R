#
# Solve a basic RBC model with BMR using gensys
#
rm(list=ls())
library(BMR)
#
beta      <- 0.99
alpha     <- .33
delta     <- .015
eta       <- 1.0
rho       <- 0.95
sigmaTech <- 1
#
RSS  <- 1/beta
YKSS <- (RSS + delta - 1)/alpha
IKSS <- delta
IYSS <- ((alpha*delta)/(RSS + delta - 1))
CYSS <- 1 - IYSS
#
Gam62 = alpha*YKSS/(RSS)
Gam63 = alpha*YKSS/(RSS)
#
#                       c_t,     y_t,     k_t,     n_t,     r_t,     i_t,     a_t,   c_t+1,   r_t+1
Gamma0 = rbind(c(         1,       0,       0,       0,       0,       0,       0,      -1,   1/eta),
               c(     -CYSS,       1,       0,       0,       0,   -IYSS,       0,       0,       0),
               c(         0,       1,       0,-1+alpha,       0,       0,      -1,       0,       0),
               c(         0,       0,       1,       0,       0,   -IKSS,       0,       0,       0),
               c(      -eta,       1,       0,      -1,       0,       0,       0,       0,       0),
               c(         0,   Gam62,       0,       0,      -1,       0,       0,       0,       0),
               c(         0,       0,       0,       0,       0,       0,       1,       0,       0),
               c(         1,       0,       0,       0,       0,       0,       0,       0,       0),
               c(         0,       0,       0,       0,       1,       0,       0,       0,       0))

#                      [c_t,     y_t,     k_t,     n_t,     r_t,     i_t,     a_t,   c_t+1,   r_t+1]
Gamma1 = rbind(c(         0,       0,       0,       0,       0,       0,       0,       0,       0),
               c(         0,       0,       0,       0,       0,       0,       0,       0,       0),
               c(         0,       0,   alpha,       0,       0,       0,       0,       0,       0),
               c(         0,       0, 1-delta,       0,       0,       0,       0,       0,       0),
               c(         0,       0,       0,       0,       0,       0,       0,       0,       0),
               c(         0,       0,   Gam63,       0,       0,       0,       0,       0,       0),
               c(         0,       0,       0,       0,       0,       0,     rho,       0,       0),
               c(         0,       0,       0,       0,       0,       0,       0,       1,       0),
               c(         0,       0,       0,       0,       0,       0,       0,       0,       1))
# 
C = matrix(0,9,1)
#
Psi = matrix(0,9,1)
Psi[7,1] = 1
#     
Pi = matrix(0,9,2)
Pi[8,1] = 1
Pi[9,2] = 1
#
# Solve first using the `gensys' function
dsgetest1 <- gensys(Gamma0,Gamma1,C,Psi,Pi)
#
# Next solve using `SDSGE'
mats <- list()
mats$Gamma0 <- Gamma0; mats$Gamma1 <- Gamma1; mats$C <- C; mats$Psi <- Psi; mats$Pi <- Pi
dsgetest2 <- SDSGE(mats)
dsgetest3 <- SDSGE(mats,type=1)
#
# all 3 dsgetest* should give the same solution, but with different 'class' labels
#
varnames=c("Consumption","Output","Capital","Labour","Interest","Investment","Technology")
dsgetestirf1 <- IRF(dsgetest1,sigmaTech,30,varnames=varnames,save=FALSE)
dsgetestirf2 <- IRF(dsgetest2,sigmaTech,30,varnames=varnames,save=FALSE)
dsgetestsim <- DSGESim(dsgetest1,1,200,200,1122,hpfiltered=FALSE)
#
#
#END