
bayes <- function(data,nox1=FALSE,...){
    sdat <- with(data, list(
                           nctl=sum(1-Z),
                           ntrt=sum(Z),
                           x1t=x1[Z==1],
                           x1c=x1[Z==0],
                           x2t=x2[Z==1],
                           x2c=x2[Z==0],
                           Ytrt=Y[Z==1],
                           Yctl=Y[Z==0],
                           St=S[Z==1]
                       )
                 )
    if(nox1){
	output <- capture.output(fit <- sampling(stanModNox1,data=sdat))#,...))
    } else{
                                       #fit1 <-
        output <- capture.output(fit <- sampling(stanPSmod,data=sdat))#,...))
    }

    summary(fit, par=c('eff0','eff1','effDiff'))$summary
}


### mu00=0
makeDat <- function(n,mu01,mu10,mu11,b1,errDist,intS,intZ,oneVarInt=TRUE,debug=FALSE){

    if(debug){
      print("debugging--parameter values set to defaults")
      if(missing(n)) n <- 1000
      if(missing(mu01)) mu01 <- 0.3
      if(missing(mu10)) mu01 <- 0.3
      if(missing(mu11)) mu11 <- 0.3
      if(missing(b1)) b1 <- 0.5
      if(missing(errDist)) errDist <- "norm"
      if(missing(intS)) intS <- FALSE
      if(missing(intZ)) intZ <- FALSE
      }

    x1 <- rnorm(2*n)
    x2 <- rnorm(2*n)
    x3 <- if(errDist=='norm'){
        rnorm(2*n)
     } else if(errDist=="lognorm") {
        exp(rnorm(2*n,0,sqrt(log((1+sqrt(5))/2))))
      } else runif(2*n,-sqrt(12)/2,sqrt(12)/2)

    x1 <- x1-mean(x1)
    x2 <- x2-mean(x2)
    x3 <- x3-mean(x3)

    psTrue <- plogis(b1*(x1+x3)-b1*x2)

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)

    error <- if(errDist=='norm') rnorm(2*n,0,sqrt(0.5))
             else if(errDist=="lognorm") exp(rnorm(2*n,0,sqrt(log((1+sqrt(3))/2))))
             else runif(2*n,-sqrt(6)/2,sqrt(6)/2)
    error <- error-mean(error)

    Yc <- sqrt(1/6)*(x1+x2+x3)+mu01*S+error
    if(intS) Yc <- Yc-sqrt(1/6)/4*(x1+x2)+S*sqrt(1/6)/2*(x1+x2)

    Yt <- Yc+mu10+(mu11-mu01-mu10)*S
    if(intZ)
	Yt <- Yt+ sqrt(1/6)/2*(if(oneVarInt) x1 else x1+x2)


    Y <- ifelse(Z==1,Yt,Yc)

    dat <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2,Strue=S)
    attr(dat,'trueEffs') <- c(
        S0=mean(Yt[S==0])-mean(Yc[S==0]),
        S1=mean(Yt[S==1])-mean(Yc[S==1])
    )
    attr(dat,'facs') <- as.data.frame(as.list(match.call())[-1])

    dat

}


simOneBayes <- function(dat,nox1=FALSE,geepers=TRUE,pmm=TRUE,psweight=TRUE){

    if(geepers) mest <- if(nox1) effs(dat,covFormU=~x2) else effs(dat) else mest <- NULL
    BAYES <- if(pmm) bayes(dat,nox1=nox1,chains=2,iter=3000,warmup=1000) else NULL
    PSW <-  if(psweight) pswSim(dat,nox1=nox1) else NULL

    list(
        true=attr(dat,'trueEffs'),
        mest=mest,
        psw=PSW,
        bayes=BAYES,
        facs=attr(dat,'facs')
    )
}


oneCase <- function(nsim,ext,ncores,cl=NULL, facs,oneVarInt=TRUE,nox1=FALSE,geepers=TRUE,pmm=TRUE,psweight=TRUE){


    if(!is.null(cl)) clusterExport(cl,list="facs",envir=environment())

     datasets <-
     if(is.null(cl)) mclapply(1:nsim,function(i) do.call("makeDat",facs),mc.cores=ncores)
     else parLapply(cl, 1:nsim,function(i) do.call("makeDat",facs))

    save(datasets,file=paste0('simData/dat',ext,'.RData'))

    estimate(datasets,nox1=nox1,geepers=geepers,pmm=pmm,psweight=psweight,cl=cl)

}

estimate <- function(datasets,nox1=FALSE,geepers=TRUE,pmm=TRUE,psweight=TRUE,cl=NULL){
	 nsim <- length(datasets)

     if(!is.null(cl)){
        clusterExport(cl,"datasets",envir=environment())

        pb <- txtProgressBar(max = nsim, style = 3)
        progress_fun <- function(nn) setTxtProgressBar(pb, nn)
        opts <- list(progress = progress_fun)
    }

    startTime <- Sys.time()
    res <-
        if(is.null(cl)){
            lapply(
                datasets,
                function(dat)
                    try(simOneBayes(dat,nox1=nox1,geepers=geepers,pmm=pmm,psweight=psweight)))

	     } else{
	     	    foreach(i = 1:nsim,
                    .packages=c("sandwich","dplyr","rstan"),
                    .options.snow = opts) %dopar% {
                        try(simOneBayes(datasets[[i]],nox1=nox1,geepers=geepers,pmm=pmm,psweight=psweight))
                    }
         }

    time <- Sys.time() - startTime

  attr(res,"time") <- time
  res
}


fullsim <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),
                    mu10=c(0,.3),
                    mu11=.3,
                    errDist=c('norm','lognorm','unif'),
                    b1s=c(0,0.3,0.5),
                    ext='',
                    intS=c(TRUE,FALSE),
                    intZ=c(TRUE,FALSE),
                    se=TRUE,
                    ncores=8,
		    cl=NULL,
		    start=1,
		    oneVarInt=TRUE,
		    nox1=FALSE,
		    geepers=TRUE,pmm=TRUE,psweight=TRUE
                    ){

    cases=expand.grid(
		n=ns,
		mu01=mu01,
		mu10=mu10,
		mu11=mu11,
		errDist=errDist,
                b1=b1s,
                intS=intS,
                intZ=intZ,
		oneVarInt=oneVarInt,
		stringsAsFactors=FALSE)

    cat(nrow(cases),' conditions\n')

    save(cases,file=paste0('simResults/cases',ext,'.RData'))

    for(i in start:nrow(cases)){
    	  cat(round(i/nrow(cases)*100),'%\n')
	  facs <- cases[i,]
    	  res <- oneCase(nsim=nsim,ext=paste0(i,ext),
                         ncores=ncores,cl=cl,facs=facs,oneVarInt=oneVarInt,nox1=nox1,geepers=geepers,pmm=pmm,psweight=psweight)
	  save(res,facs,file=paste0('simResults/sim',i,ext,'.RData'))
    }

    return(0)
}



pswSim=function(dat,psMod,nox1){

  if(missing(psMod)){
	if(nox1){
		psMod=glm(S~x2,data=dat,subset=Z==1,family=binomial)
		} else psMod=glm(S~x1+x2,data=dat,subset=Z==1,family=binomial)
  }
  dat0=subset(dat,Z==0)
  dat1=subset(dat,Z==1)
  dat0$ps=predict(psMod,dat0,type='response')


  muc0=with(dat0,sum(Y*(1-ps))/sum(1-ps))
  muc1=with(dat0,sum(Y*ps)/sum(ps))

  mut0=with(dat1,mean(Y[S==0]))
  mut1=with(dat1,mean(Y[S==1]))

  c(eff0=mut0-muc0,
    eff1=mut1-muc1,
    attributes(dat)$trueEffs
    )
}

