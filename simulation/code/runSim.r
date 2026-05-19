library(rstan)
library(parallel)
library(sandwich)
library(snow)
library(doSNOW)
library(utils)
library(tidyverse)
library(broom)
library(kableExtra)
library(gridExtra)
library(ggpubr)
library(tikzDevice)



rstan_options(auto_write = TRUE)

##########################################################################################
########### run the simulation
##########################################################################################



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

#### for parallel computing
if(ncore>1){
    print("making cluster")
    cl <- makeCluster(ncore)
    registerDoSNOW(cl)
    clusterEvalQ(cl,library(dplyr))
    clusterEvalQ(cl,library(rstan))
    clusterEvalQ(cl,rstan_options(auto_write=TRUE))
    clusterEvalQ(cl, stanModNox1 <- stan_model('code/psNoX1.stan',model_name="stanModNox1"))
    clusterEvalQ(cl, stanPSmod <- stan_model("code/ps.stan",model_name="stanPSmod"))
    clusterEvalQ(cl,source("../regression.r"))
    clusterEvalQ(cl,source('code/simFuncs.r'))
} else{
    stanModNox1 <- stan_model('code/psNoX1.stan',model_name="stanModNox1")
    stanPSmod <- stan_model("code/ps.stan",model_name="stanPSmod")
    cl <- NULL
}


source('code/simFuncs.r')

source('../regression.r')

set.seed(613)
fullsim(nrep,ns=c(500),mu01=c(0),mu10=0,b1s=c(0,0.3,0.5),
	ext='',cl=cl,ncores=ncore,start=1)



################################
### look at wider ranges of n and alpha

fullsim(nrep,ns=c(100,200,300,400,600,700,800,900,1000),errDist='norm',b1s=.5,intS=FALSE,intZ=FALSE,mu01=0,mu10=0,ext='ns',cl=cl,
        ncores=ncore)

fullsim(nrep,ns=500,errDist='norm',b1s=c(.1,.2,.4,seq(.6,1,.1)),intS=FALSE,intZ=FALSE,
  mu01=0,mu10=0,ext='b1s',cl=cl, ncores=ncore)


################################
### try estimation excluding x1
### use the same simulated datasets as before

data=list.files("simData")
data=grep("dat[1-9]*\\.RData",data,value=TRUE)

for(dd in data){
       print(dd)
       load(paste0("simData/",dd))
       res=estimate(datasets,nox1=TRUE,pmm=TRUE,psweight=TRUE,cl=cl)
       save(res,file=paste0("simResults/",gsub("dat","sim",gsub("\\.RData","nox1.RData",dd))))
}


stopCluster(cl)




##########################################################################################
########### read in the results, pre-process, and make tables and figures
##########################################################################################


 source('code/readSim.r')

 source('code/displaySimPublish.r')

#if(!is.null(cl))

