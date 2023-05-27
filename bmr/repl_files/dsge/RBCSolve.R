#
# Solve a basic RBC model with BMR using Uhlig's method
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
A <- matrix(c( 0,
               0,
               1,
               0,
               0 ))

#
B <- matrix(c(                0,
                         -alpha,
                     -(1-delta),
                              0,
               (alpha/RSS)*YKSS))

#              consumption,              output,          labor,    interest,   investment
C <- rbind(c(        -CYSS,                   1,              0,           0,      -IYSS), 
           c(            0,                   1,     -(1-alpha),           0,          0),
           c(            0,                   0,              0,           0,      -IKSS),      
           c(          eta,                  -1,              1,           0,          0),
           c(            0,   -(alpha/RSS)*YKSS,              0,           1,          0))
D <- matrix(c(  0,
               -1,
                0,
                0,
                0))
F <-  matrix(c(0))
G <- matrix(c(0))
H <- matrix(c(0))
J <- t(matrix(c( -1,  0,  0,  (1/eta),  0)))
K <- t(matrix(c( 1,   0,  0,  0,  0)))
L <- matrix(c(0)) 
M <- matrix(c(0)) 
N <- matrix(rho)
Sigma <- matrix(sigmaTech)
#
# Solve first using the `uhlig' function
dsgetest1 <- uhlig(A,B,C,D,F,G,H,J,K,L,M,N)
#
# Next solve using `SDSGE'
mats <- list()
mats$A <- A; mats$B <- B; mats$C <- C; mats$D <- D
mats$F <- F; mats$G <- G; mats$H <- H
mats$J <- J; mats$K <- K; mats$L <- L; mats$M <- M; mats$N <- N
dsgetest2 <- SDSGE(mats)
dsgetest3 <- SDSGE(mats,type=2)
#
# all 3 dsgetest* should give the same solution, but with different 'class' labels
#
varnames=c("Capital","Consumption","Output","Labour","Interest","Investment","Technology")
dsgetestirf1 <- IRF(dsgetest1,sigmaTech,30,varnames=varnames,save=FALSE)
dsgetestirf2 <- IRF(dsgetest2,sigmaTech,30,varnames=varnames,save=FALSE)
dsgetestsim <- DSGESim(dsgetest1,1,200,200,1122,hpfiltered=FALSE)
#
#
#END