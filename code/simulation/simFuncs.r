#stanMod <- rstan::stan_model('ps.stan')#,auto_write=TRUE)

bayes <- function(data,returnFit=FALSE,...){
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
                                        #fit1 <-
    output <- capture.output(fit <- stan('code/ps.stan',data=sdat,...))#,...))
    if(returnFit) return(fit)
    summary(fit, par=c('eff0','eff1','effDiff'))$summary
}

pstrata <- function(data,...){

    psobj <- PSObject(
        S.formula=Z+S~x1+x2,
        Y.formula=Y~x1+x2,
        Y.family=gaussian(link="identity"),
        data=data,
        strata = c(nt = "00", co = "01"),
        ER=c(nt=FALSE,co=FALSE)
    )
    standata <- make_standata(psobj)
    output <-
        capture.output(
            stansamp <- post_samples <- PSSample(
                    "code/pstrata.stan",
                    data = standata))
    ps1 <-    list(
        PSobject = psobj,
        post_samples = post_samples
    )
    class(ps1) <- "PStrata"
    ps1%>%PSOutcome()%>%PSContrast(Z=TRUE)%>%summary("matrix")
}


makeDatPois <- function(n,mu01,mu10,mu11,b1,intS,intZ,oppSignsPS=TRUE,
                        intZfunc=\(x1,x2) sqrt(1/6)/2*(x1+x2)){

    x1=exp(rnorm(2*n,0,sqrt(s2)))
    x2=exp(rnorm(2*n,0,sqrt(s2)))
    x3=exp(rnorm(2*n,0,sqrt(s2)))

    x1 <- runif(n,0,1/4)
    x2 <- runif(n,0,1/4)
    x3 <- runif(n,0,1/4)


    linPred <- if(oppSignsPS) b1*(x1+x3)-b1*x2 else b1*(x1+x3+x2)
    fout <- rep(c(0,1),n)
    intercept <- coef(glm(fout~1,family=binomial,offset=b1*(x1+x2+x3)/3))[1]
    psTrue <- plogis(b1*(x1+x2+x3)/3+intercept)

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)

    #beta <- 1/sqrt(6)#1/(6*exp(s2/2))

    beta <- mean(x1+x2+x3)/var(x1+x2+x3)

    lambdaC <- beta*(x1+x2+x3)+mu01*S
    if(intS) lambdaC <- lambdaC-beta/4*(x1+x2)+S*beta/2*(x1+x2)

    lambdaT <- lambdaC+mu10+(mu11-mu01-mu10)*S
    if(intZ) lambdaT <- lambdaT+intZfunc(x1,x2)

    Yt <- rpois(2*n,lambdaT)
    Yc <- rpois(2*n,lambdaC)
    Y <- Yc*(1-Z)+Yt*Z
    sig <- sd(Y)

    dat <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2,Strue=S)
    attr(dat,'trueEffs') <- c(
        S0=mean(Yt[S==0]-Yc[S==0]),
        S1=mean(Yt[S==1]-Yc[S==1])
    )#/sig
    facs <- as.list(match.call())[-1]
    facs$intZfunc <- NULL
    attr(dat,'facs') <- facs
    attr(dat, "intZfunc") <- intZfunc

    dat
}


s <- log(37)
m <- -s/2-log(6)

### mu00=0
makeDat <- function(n,mu01,mu10,mu11,b1,errDist,intS,intZ,oppSignsPS=TRUE,
                    intZfunc=\(x1,x2) sqrt(1/6)/2*(x1+x2)){

    if(errDist=="pois")
        return(makeDatPois(n=n,mu01=mu01,mu10=mu10,mu11=mu11,b1=b1,intS=intS,intZ=intZ,
                           oppSignsPS=oppSignsPS,intZfunc=intZfunc))

    x1 <- rnorm(2*n)
    x2 <- rnorm(2*n)
    x3 <- if(is.character(errDist)){
              if(errDist=='norm'){
                  rnorm(2*n)
              } else if(errDist=="lognorm") {
                  exp(rnorm(2*n,0,sqrt(log((1+sqrt(5))/2))))
              } else if(errDist=="unif"){
                  runif(2*n,-sqrt(12)/2,sqrt(12)/2)
              }
          } else errDist(2*n)

    x1 <- x1-mean(x1)
    x2 <- x2-mean(x2)
    x3 <- (x3-mean(x3))/sd(x3)

    psTrue <- if(oppSignsPS) plogis(b1*(x1+x3)-b1*x2) else plogis(b1*(x1+x3+x2))

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)

    error <- if(is.character(errDist)){
                 if(errDist=='norm') rnorm(2*n,0,sqrt(0.5))
                 else if(errDist=='mix') c(rnorm(3*n/2,-1/3,sqrt(1/6)),rnorm(n/2,1,sqrt(1/6)))
                 else if(errDist=="lognorm") exp(rnorm(2*n,0,sqrt(log((1+sqrt(3))/2))))
                 else if(errDist=="unif") runif(2*n,-sqrt(6)/2,sqrt(6)/2)
             }  else errDist(2*n)
    error <- (error-mean(error))/sd(error)*sqrt(0.5)

    Yc <- sqrt(1/6)*(x1+x2+x3)+mu01*S+error
    if(intS) Yc <- Yc-sqrt(1/6)/4*(x1+x2)+S*sqrt(1/6)/2*(x1+x2)

    #x0 <- mean(x1[S==0]+x2[S==0])

    Yt <- Yc+mu10+(mu11-mu01-mu10)*S
    if(intZ) Yt <- Yt+intZfunc(x1,x2)

    Y <- ifelse(Z==1,Yt,Yc)

    ## if(!norm){
    ##     Y <- round(Y)
    ##     #Y[Y< -4] <- 4
    ##     #Y[Y> 4] <- 4
    ## }

    dat <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2,Strue=S)
    attr(dat,'trueEffs') <- c(
        S0=mean(Yt[S==0])-mean(Yc[S==0]),
        S1=mean(Yt[S==1])-mean(Yc[S==1])
    )
    facs <- as.list(match.call())[-1]
    facs$intZfunc <- NULL
    attr(dat,'facs') <- facs
    attr(dat, "intZfunc") <- intZfunc

    dat

}


simOneBayes <- function(dat){

    mest <- effs(dat)
    BAYES <- bayes(dat,chains=2,iter=3000,warmup=1000)
    PSW <- psw(dat)

    list(
        true=attr(dat,'trueEffs'),
        mest=mest,
        psw=PSW,
        bayes=BAYES,
        facs=attr(dat,'facs')
    )
}


simOne <- function(dat,estimators=list(mest=effs,bayes=bayes,psw=psw)){

    out <- lapply(estimators, \(fun) fun(dat))
    out$true <- attr(dat,'trueEffs')
    out$facs <- attr(dat,"facs")

    out
}
  ##   list(

##         mest=mest,
##         psw=PSW,
##         bayes=BAYES,
##         facs=attr(dat,'facs')
##     )
## }


oneCase <- function(nsim,ext,ncores,cl=NULL, facs,estimators,intZfunc){ #n,mu00,mu01,mu10,mu11,gumb,b1,cl){

#    print(Sys.time())

    if(!is.null(cl)) clusterExport(cl,list="facs",envir=environment())

    facs=c(facs,intZfunc=intZfunc)

     datasets <-
     if(is.null(cl)) mclapply(1:nsim,function(i) do.call("makeDat",facs),mc.cores=ncores)
     else parLapply(cl, 1:nsim,function(i) do.call("makeDat",facs))

    save(datasets,file=paste0('simData/dat',ext,'.RData'))

    if(!is.null(cl)){
        clusterExport(cl,"datasets",envir=environment())
        clusterExport(cl,"estimators",envir=environment())

        pb <- txtProgressBar(max = nsim, style = 3)
        progress_fun <- function(nn) setTxtProgressBar(pb, nn)
        opts <- list(progress = progress_fun)
    }

    startTime <- Sys.time()
    res <-
        if(is.null(cl)){
	    	mclapply(
		    	datasets,
                	function(dat) try(simOne(dat,estimators=estimators)),
                	mc.cores=ncores
            		)
	     } else{
	     	    foreach(i = 1:nsim,
                    .packages=c("sandwich","dplyr","rstan"),
                    .options.snow = opts) %dopar% {
                        try(simOne(datasets[[i]],estimators=estimators))
                    }
         }
#                  parLapply(cl,datasets, function(dat) try(simOneBayes(dat)))

    time <- Sys.time() - startTime
#    print(time)

  attr(res,"time") <- time
  res
}

fullsim <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),#sepTs=c(TRUE,FALSE),
                    mu10=c(0,.3),#sepCs=c(TRUE,FALSE),
                    mu11=.3,#effs=c(TRUE,FALSE),
                    errDist=c('norm','lognorm','unif'),
                    b1s=c(0,0.3,0.5),
                    oppSignsPS=TRUE,
                    ext='',
                    intS=c(TRUE,FALSE),
                    intZ=c(TRUE,FALSE),
                    ncores=8,
                    cl=NULL,
                    start=1,
                    estimators=list(mest=effs,bayes=bayes,psw=psw),
                    intZfunc=\(x1,x2) sqrt(1/6)/2*(x1+x2)
                    ){

    cases=expand.grid(
		n=ns,
		mu01=mu01,
		mu10=mu10,
		mu11=mu11,
    oppSignsPS=oppSignsPS,
		errDist=errDist,
    b1=b1s,
    intS=intS,
    intZ=intZ,
    #intZfunc=intZfunc,
		stringsAsFactors=FALSE)

    cat(nrow(cases),' conditions\n')

    save(cases,file=paste0('simResults/cases',ext,'.RData'))

    for(i in start:nrow(cases)){
    	  cat(round(i/nrow(cases)*100),'%\n')
	  facs <- cases[i,]
    	  res <- oneCase(nsim=nsim,ext=paste0(i,ext),
                         ncores=ncores,cl=cl,facs=facs,estimators=estimators,intZfunc)
	  save(res,facs,file=paste0('simResults/sim',i,ext,'.RData'))
    }

    return(0)
}

fullsimJustM <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),#sepTs=c(TRUE,FALSE),
                    mu10=c(0,.3),#sepCs=c(TRUE,FALSE),
                    mu11=.3,#effs=c(TRUE,FALSE),
                    gumbs=c(TRUE,FALSE),
                    b1s=c(0,0.2,0.5,1),
                    ext='',
                    se=TRUE,
                    cl=NULL,
                    Bayes=FALSE,
		    start=1
                    ){

    cases=expand.grid(ns,mu01,mu10,mu11,gumbs,b1s)
    names(cases) <- c('n','mu01','mu10','mu11','gumb','b1')

    if(nsim==0) return(cases)

    cat('% done:')
    list(cases=cases,
         res=lapply(
             1:nrow(cases),
             function(i){
                 cat(round(i/nrow(cases)*100))
                 replicate(nsim,effs(do.call("makeData",cases[i,])))
             }
         )
         )
}


psw=function(dat,psMod){

  if(missing(psMod)) psMod=glm(S~x1+x2,data=dat,subset=Z==1,family=binomial)

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


xSim <- function(){
     n <- 10000
 mu01=0.2
 mu10=0
 mu11=0.5
 b1=1
 errDist='norm'
 x1 <- rnorm(2*n)
 x2 <- rnorm(2*n)
 x3 <- if(errDist=='norm') rnorm(2*n) else runif(2*n,-sqrt(12)/2,sqrt(12)/2)
 x1 <- x1-mean(x1)
 x2 <- x2-mean(x2)
 x3 <- x3-mean(x3)
 psTrue <- plogis(b1*(x1+x3)-b1*x2)
 S <- rbinom(2*n,1,psTrue)
 Z <- rep(c(1,0),n)
 x12 <- x1+x2
     c(mean(x12),
       mean(x12[S==1]),
       mean(x12[S==0]))
}


