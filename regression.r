library(sandwich)

DBN <- #"none"
## "n1"
"n"


### I want a bunch of attributes for the VCV matrix, but I don't
### want them to print every time!
print.vcv <- function(x,...) print(x[1:NROW(x),1:NCOL(x)],...)

print.geepers <- function(x,...) print(effsFromFit(x),...)

### convenience function for functions that return lists
### why cant R be more like python in this regard??
Attach <- function(.list,.names=names(.list),env=environment()){
    if(is.null(.names))
        .names <- paste0(
            as.character(
                as.list(match.call())$.list)[1],
            seq(length(.list)))

    for(i in 1:length(.list)){
        assign(.names[i],.list[[i]], envir=parent.frame())
    }

}

auc <- function(x,y)
  if(length(unique(y))==2 & length(y)==length(x)){
    wilcox.test(x~y)$statistic/(sum(y)*(sum(1-y)))
  }else
    wilcox.test(x,y)$statistic/(legnth(x)*length(y))

aucMod <- function(mod)
  auc(mod$linear,1-mod$y)



AUC1 = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(auc)
}
AUCmod <- function(mod){
    auc(mod$linear,1-mod$y)
}

datNames <- function(data,trt='Z',out='Y',use='S',clust=NULL,block=NULL){
  data$Y <- as.numeric(data[[out]])
  data$Z <- as.numeric(data[[trt]])
  data$S <- as.numeric(data[[use]])
  data$clust <- if(is.character(clust)) data[[clust]] else clust

  data$S[data$Z==0] <- NA

  if(!is.null(block)){
    data$block <- as.factor(data[[block]])
    data <- data%>%
      group_by(block)%>%
      mutate(
        Zadj=Z-mean(Z),
        Yadj=Y-mean(Y)
      )
  } else{
    data <- mutate(data,
                   Z=Z,
                   Yadj=Y)
  }
  data
}

getPointEst <- function(outMod,covFormY,data){

    xbar0 <- colMeans(model.matrix(covFormY,data=data[data$Z==1&data$S==0,])[,-1,drop=FALSE])
    xbar1 <- colMeans(model.matrix(covFormY,data=data[data$Z==1&data$S==1,])[,-1,drop=FALSE])

    stopifnot(all.equal(names(xbar0),names(xbar1)))

    coefs <- coef(outMod)

    beta2 <- coefs["Z"]
    beta3 <- coefs["Z:Sp"]

    gamma3 <- coefs[paste0("Z:",names(xbar0))]
    gamma4 <- coefs[paste0("Z:Sp:",names(xbar0))]

    if(is.na(gamma3[1])) gamma3 <- 0
    if(is.na(gamma4[1])) gamma4 <- 0

    out <- list(
        eff0=unname(beta2+sum(gamma3*xbar0)),
        eff1=unname(beta2+beta3+sum((gamma3+gamma4)*xbar1))
    )

    out$diff <- out$eff1-out$eff0
    out

}

### computes/returns regression models
pointEst <- function(data,covFormU=~x1+x2,covFormY=covFormU,intSx=NULL,psMod=NULL,bayes=FALSE){

    if(bayes) stopifnot(require(arm))

  data$S[data$Z==0] <- NA

    if(is.null(psMod)){
        psMod <- if(bayes){
                     bayesglm(update(covFormU,S~.),data=subset(data,!is.na(S)&Z==1),family=binomial)
                 } else glm(update(covFormU,S~.),data=subset(data,!is.na(S)),family=binomial)
    } else covFormU <- formula(psMod)[c(1,3)]

    if(is.null(covFormY)) covFormY <- covFormU


    attr(psMod,'auc') <- AUCmod(psMod)

    ps <- predict(psMod,data,type='response')
    if('Sp'%in%names(data)) warning('replacing Sp')
    data <- within(data,Sp <- ifelse(Z==1&!is.na(S),S,ps))
  outForm <- update(covFormY,Y~Z*Sp+.)

#  outForm <- update(covFormY,Y~Z*Sp+.)

  if(!is.null(intSx))
    outForm <- update(outForm,as.formula(paste0('.~.+Sp:(',paste(intSx,collapse='+'),')')))


  if('block'%in%names(data)) if(!is.null(data$block))
      outForm <- update(outForm,.~.+block)

    outMod <- if(bayes){
                  bayesglm(update(covFormY,Y~Z*Sp*.),data=data,prior.scale=1,prior.df=Inf,scaled=FALSE)
              } else lm(outForm,data=data)

    estimates <- getPointEst(outMod,covFormY,data)

    list(psMod=psMod,outMod=outMod,estimates=estimates)
}

intEst0 <- function(dat){
  Attach(pointEst(dat,int=TRUE))

  coef(outMod)['Z']+
    crossprod(coef(outMod)[c('Z:x1','Z:x2')],colMeans(subset(dat,Z==1&S==0,select=c('x1','x2'))))
}

### estimating equations for stage 1 model (PS)
EFps <- function(psMod,data){
  efPStmp <- estfun(psMod)
  efPS <- matrix(0,nrow(data),ncol(efPStmp))
  if(is.null(rownames(data))){
    efPS[as.numeric(rownames(efPStmp)),] <- efPStmp
  } else efPS[match(rownames(efPStmp),rownames(data)),] <- efPStmp
    efPS
}

### estimates bread and meat for sandwich variance
sandwichMats <- function(psMod,outMod,data,clust=NULL,int=any(grepl(":x",names(coef(outMod))))){
    #data <- model.frame(outMod)
    ### estimating equations
    efPS <- EFps(psMod,data) ## psi1
    efOut <- estfun(outMod) ## psi2

    out <- list(
        a11inv = bread(psMod)/(if(DBN=="none") nobs(psMod) else 1),
        a22inv = bread(outMod)/(if(DBN=="none") nobs(outMod) else 1),
        a21 = A21(psMod,outMod,data)/(if(DBN=="none") 1 else if(DBN=="n") nrow(data) else nobs(psMod)),
        b11=(if(is.null(clust)) meat(psMod) else meatCL(psMod,cluster=data[rownames(model.frame(psMod)),"clust"]))*
            (if(DBN=="none") nobs(psMod) else 1),
        b22 = (if(is.null(clust)) meat(outMod) else meatCL(outMod,cluster=data$clust))*
            (if(DBN=="none") nobs(outMod) else 1)
    )
    out <- within( out,b12 <-
        if(int){
            matrix(0,nrow(b11),ncol(b22))
        } else if(is.null(clust)) crossprod(efPS,efOut)/
                                  (if(DBN=="none") 1 else if(DBN=="n") nrow(data) else nobs(psMod))
        else crossprod(apply(efPS,2L,rowsum,data$clust),
                       apply(efOut,2L,rowsum,data$clust))/
             (if(DBN=="none") 1 else if(DBN=="n") nrow(data) else nobs(psMod))
        )
    out
}

### combining smaller matrices into bigger one
bigM <- function(m11,m12,m22)
    rbind(
        cbind(m11,m12),
        cbind(t(m12),m22)
    )

### computes sandwich vcov for regressions
vcvPS <- function(psMod,outMod,data,clust=NULL,int=any(grepl(":x",names(coef(outMod))))){
    Attach(
        sandwichMats(
            psMod=psMod,
            outMod=outMod,
            data=data,
            clust=clust,
            int=int))

    A <- rbind(
        cbind(solve(a11inv),matrix(0,nrow(a11inv),ncol(a22inv))),
        cbind(a21,solve(a22inv))
    )
    colnames(A) <- rownames(A)
    ##bigM(solve(a11inv),t(a21),solve(a22inv))

    B <- bigM(b11,b12,b22)

    n <- nrow(model.frame(outMod))

    vcvFull <- solve(A)%*%B%*%t(solve(A))/(if(DBN=="none") 1 else nrow(data))
    main <- a22inv%*%b22%*%t(a22inv)
    vcv <- vcvFull[-(1:nrow(b11)),-(1:nrow(b11))]
    dimnames(vcv) <- dimnames(main)
    correction <- vcv-main

    attr(vcv,"full") <- vcvFull
    attr(vcv,"bread") <- A
    attr(vcv,"meat") <- B
    attr(vcv,"main") <- main
    attr(vcv,"correction") <- correction

    class(vcv) <- "vcv"

    vcv
}



### wrapper function for estimating regressions+vcov
est <- function(data,covFormU=~x1+x2,covFormY=NULL,psMod=NULL,clust=NULL,intSx=NULL,
                trt='Z',out='Y',use='S',block=NULL,bayes=FALSE){

  #if(!missing(psMod))

  data <- datNames(data,trt=trt,out=out,use=use,block=block,clust=clust)

  Attach(pointEst(data=data,covFormU=covFormU,covFormY=covFormY,intSx=intSx,psMod=psMod,bayes=bayes))

  vcv <- vcvPS(psMod,outMod,data=data,clust=clust)

  out <- list(outMod=outMod,psMod=psMod,vcv=vcv,estimates=estimates)
  class(out) <- c('geepers',class(out))
  out
}

### estimates effects of interest, starting from est() output
effsFromFit <- function(ests){
    estimates <- ests$estimates

    vcv <- ests$vcv
    ddd <- diag(vcv)
    vars <- list(
        eff0 <- ddd['Z'],
        eff1 <- ddd['Z']+ddd['Z:Sp']+2*vcv['Z','Z:Sp'],
        diff <- ddd['Z:Sp']
    )
    out <-
        cbind(
            estimates=unlist(estimates),
            SE=sqrt(unlist(vars))
        )
    attr(out,'auc') <- attr(ests$psMod,'auc')
    out
}

### estimates effects of interest, starting from data
effs <- function(data,covFormU=~x1+x2,covFormY=covFormU,intSx=NULL,
                 trt='Z',out='Y',use='S',block=NULL){
    ests <- est(data=data,covFormU=covFormU,covFormY=covFormY,
                intSx=intSx,trt=trt,out=out,use=use,block=block)

    effsFromFit(ests)
}

### estimates lower left of A (bread) matrix in A^{-1}BA^{-t}
### INTERACTION W SP ONLY WORKS IF NO MISSING S IN TRT GROUP
A21 <- function(psMod,outMod,data){

  intSx <- any(grepl('Sp:',names(coef(outMod))))

  if(intSx&any(is.na(data$S[data$Z==1])))
    warning("INTERACTION W SP ONLY WORKS IF NO MISSING S IN TRT GROUP")

  risp <- data$Z==0|is.na(data$S)#,1,0)
  Z <- data$Z[risp]
  Y0 <- data$Y[risp]
  aW <- predict(psMod,data,type='link')[risp]  ### X'theta
  p <- family(psMod)$linkinv(aW)
  q <- family(psMod)$mu.eta(aW)      ## dp/d(aW) so that dp/dalpha=qW'

  W0 <- model.matrix(update(formula(psMod),Y~.),data=data)[risp,]

  X0 <- model.matrix(outMod)[risp,
                               -c(which(names(coef(outMod))%in%c('Z','Sp',"Z:Sp")),
                                  #grep('Z\\:|Sp\\:|\\:Z|\\:Sp',names(coef(outMod))),
                                  which(names(coef(outMod))=="(Intercept)")),drop=FALSE]

  if(intSx) V0 <- cbind(model.matrix(outMod)[risp,grep('Sp\\:',names(coef(outMod)))])

  Q <- rbind(coef(outMod)[c('(Intercept)','Z',colnames(X0))])%*%t(cbind(1,Z,X0))
  U <- rbind(coef(outMod)[c('Sp','Z:Sp')])%*%t(cbind(1,Z))
  if(intSx) U <- U+(rbind(coef(outMod)[colnames(V0)])%*%t(V0))

  AA <- rbind(
    `(Intercept)`=(-U)[1,],
    Sp=(Y0-(Q+2*p*U))[1,],
    Z=(-Z*U)[1,],
    `Z:Sp`=(Z*Y0-(Z*Q+2*p*Z*U))[1,],
    -sweep(t(X0),2,U,"*")
  )

  if(intSx) AA <- rbind(AA,(Y0-Q-2*p*U)[rep(1,ncol(V0)),]*t(V0))

  AA <- AA[names(coef(outMod)),] ## make sure the rows are in the same order as other matrices

  DD <- W0*q

  -AA%*%DD#/nrow(X0)
}


################ PSW functions
#### PSW
psw1 <- function(dat,psMod){

  dat0=subset(dat,Z==0)
  dat1=subset(dat,Z==1)
  dat0$ps=predict(psMod,dat0,type='response')


  muc0=with(dat0,sum(Y*(1-ps))/sum(1-ps))
  muc1=with(dat0,sum(Y*ps)/sum(ps))

  mut0=with(dat1,mean(Y[S==0]))
  mut1=with(dat1,mean(Y[S==1]))

  c(eff0=mut0-muc0,
    eff1=mut1-muc1,
    diff=mut1-muc1-mut0+muc0
    )
}

print.psw <- function(x, ...) print(x$coef,...)

bsSchool <- function(dat){
  datS <- split(dat,dat$SchIDPre)
  do.call("rbind",
          lapply(datS,function(x) x[sample(1:nrow(x),nrow(x),replace=TRUE),]))
}

bsInd=function(data)
  data[sample(1:nrow(data),nrow(data),replace=TRUE),]

psw <- function(dat,psMod,B=5000,verbose=TRUE,bsFun=bsInd){
  #dat <- getDat(alt)

  est=psw1(dat,psMod)

  if(verbose) step <- if(B>10) round(B/10) else 1
  bs <- matrix(nrow=B,ncol=3)
  for(i in 1:B){
    if(verbose) if(i%%step==0) cat(round(i/B*100),' ')
    datStar <- bsFun(dat)
    psModStar <- glm(formula(psMod),family=family(psMod),
                     data=subset(datStar,Z==1))
    bs[i,] <- psw1(datStar,psModStar)
  }
  if(verbose) cat('\n')
  out <-
    list(
      coef=cbind(est=est,se=apply(bs,2,sd)),
      bs=bs)
  class(out) <- 'psw'
  out
}



estfun.lmridge <- function(x,data,...){
    xmat <- model.matrix(formula(x),data=data)
    #xmat <- naresid(x$na.action, xmat)
    Kstar <- which(x$K==min(x$K[x$K>0]))
    if (any(alias <- is.na(coef(x))))
        xmat <- xmat[, !alias, drop = FALSE]
    wts <- weights(x)
    if (is.null(wts))
        wts <- 1
    res <- residuals(x)[,Kstar]
    rval <- as.vector(res) * wts * xmat+2*x$K[Kstar]*coef(x)
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    if (is.zoo(res))
        rval <- zoo(rval, index(res), attr(res, "frequency"))
    if (is.ts(res))
        rval <- ts(rval, start = start(res), frequency = frequency(res))
    return(rval)
}
