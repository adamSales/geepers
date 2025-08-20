##########################################################################################
########### run the simulation
##########################################################################################


library(dplyr)
library(rstan)
library(parallel)
library(sandwich)
library(purrr)

rstan_options(auto_write = TRUE)

#stanPSmod=stan_model("code/ps.stan",model_name="stanPSmod")


### number of repetitions
nrep <- 5000

### For running with multiple cores:

ncore <-  10

# Vector of required directories
dirs <- c("simFigs", "simResults","simData")

# Loop through and create if missing
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d)
  } 
}

#if(.Platform$OS.type=='windows' & ncore>1){
print("making cluster")
cl <- makeCluster(ncore)
clusterEvalQ(cl,library(dplyr))
clusterEvalQ(cl,library(rstan))
clusterEvalQ(cl,rstan_options(auto_write=TRUE))
clusterEvalQ(cl, stanPSmod <- stan_model("code/ps.stan",model_name="stanPSmod"))
clusterEvalQ(cl,source("code/regression.r"))
clusterEvalQ(cl,source('code/simulation/simFuncs.r'))
#} else cl <- NULL


source('code/simulation/simFuncs.r')

source('code/regression.r')


fullsim(nrep,ns=c(500),mu01=c(0),mu10=0,b1s=c(0,0.3,0.5),
	ext='',cl=cl,ncores=ncore,start=1)



################################
### look at wider ranges of n and alpha

fullsim(nrep,ns=c(100,200,300,400,600,700,800,900,1000),errDist='norm',b1s=.5,intS=FALSE,intZ=FALSE,mu01=0,mu10=0,ext='ns',cl=cl,
        ncores=ncore)

fullsim(nrep,ns=500,errDist='norm',b1s=c(.1,.2,.4,seq(.6,1,.1)),intS=FALSE,intZ=FALSE,
  mu01=0,mu10=0,ext='b1s',cl=cl, ncores=ncore)

stopCluster(cl)


##########################################################################################
########### read in the results, pre-process, and make tables and figures
##########################################################################################


# #library(tidyverse)
# library(purrr)
# library(ggplot2)
# library(broom)
# library(kableExtra)
# library(gridExtra)
# library(ggpubr)

# source('readSim.r')

# source('displaySim.r')

#if(!is.null(cl)) 

#source("code/simulation/readSim.r")
#source("code/simulation/displaySimPublish.r")
