##########################################################################################
########### run the simulation
##########################################################################################


library(dplyr)
library(rstan)
library(parallel)
library(sandwich)
library(purrr)
library(snow)
library(doSNOW)
library(utils)

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
registerDoSNOW(cl)

clusterEvalQ(cl,library(dplyr))
clusterEvalQ(cl,library(rstan))
clusterEvalQ(cl,rstan_options(auto_write=TRUE))
clusterEvalQ(cl, stanPSmod <- stan_model("code/ps.stan",model_name="stanPSmod"))
clusterEvalQ(cl, stanModNox1 <- stan_model("code/psNoX1.stan",model_name="stanModNox1"))
clusterEvalQ(cl,source("code/regression.r"))
clusterEvalQ(cl,source('code/simulation/simFuncs.r'))
#} else cl <- NULL


source('code/simulation/simFuncs.r')

source('code/regression.r')



estimate <- function(datasets,nox1=FALSE,geepers=TRUE,pmm=TRUE,psweight=TRUE){

	 nsim <- length(datasets)
#    print(Sys.time())
    
    clusterExport(cl,"datasets",envir=environment())

    pb <- txtProgressBar(max = nsim, style = 3)
    progress_fun <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress_fun)


    startTime <- Sys.time()
    res <-
        if(is.null(cl)){
	    	mclapply(
		    	datasets,
                	function(dat) try(simOneBayes(dat,nox1=nox1,geepers=geepers,pmm=pmm,psweight=psweight)),
                	mc.cores=ncores
            		)
	     } else{
	     	    foreach(i = 1:nsim,
                    .packages=c("sandwich","dplyr","rstan"),
                    .options.snow = opts) %dopar% {
                        try(simOneBayes(datasets[[i]],nox1=nox1,geepers=geepers,pmm=pmm,psweight=psweight))
                    }
         }
#                  parLapply(cl,datasets, function(dat) try(simOneBayes(dat)))

    time <- Sys.time() - startTime
#    print(time)

  attr(res,"time") <- time
  res
}


# data=list.files("simData")
# data=grep("dat[1-9]*oneVarInt\\.RData",data,value=TRUE)

# for(dd in data){
#        print(dd)
#        load(paste0("simData/",dd))
#        res=estimate(datasets,nox1=TRUE,pmm=TRUE,psweight=TRUE)
#        save(res,file=paste0("simResults/",gsub("dat","sim",gsub("\\.RData","nox1.RData",dd))))
# }

#data=grep("dat[1-9]*oneVarInt\\.RData",data,value=TRUE)
dd = "dat10oneVarInt.RData"
load(paste0("simData/",dd))
res=estimate(datasets,nox1=TRUE,pmm=TRUE,psweight=TRUE)
save(res,file="simResults/sim10oneVarIntnox1.RData")