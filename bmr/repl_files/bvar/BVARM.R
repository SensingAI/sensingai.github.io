#
# BVAR with Minnesota prior
# 
rm(list=ls())
library(BMR)
#
data(BMRVARData)
USMacroData<-USMacroData[,2:4]
#
prior <- c(0.9,0.9,0.9)
#
# Different p
#
testbvarm <- BVARM(USMacroData,prior,p=1,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
testbvarm <- BVARM(USMacroData,prior,p=2,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
testbvarm <- BVARM(USMacroData,prior,p=3,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
testbvarm <- BVARM(USMacroData,prior,p=4,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
#
#END