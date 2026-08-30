bias <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    colMeans(cbind(sss1[,3:6]-sss1[,c(1,2,1,2)],sss2))
}

rmse <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    c(
        sqrt(colMeans((sss1[,3:6]-sss1[,c(1,2,1,2)])^2)),
        colMeans(sss2)
    )
}


summ <- function(sss){
    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]
    rbind(colMeans(sss1),
          apply(sss1,2,sd),
          sqrt(row1Means((sss-sss[rep(1:2,3),])^2))
          )
    }

muEff <- function(facs,eff)
    with(as.list(facs),
         ifelse(eff==0,mu10,
         ifelse(eff==1,mu11-mu01,mu11-mu01-mu10))
         )

psw1Proc=function(psw){
  psw['effDiff']=psw['eff1']-psw['eff0']
  cbind(
    estimates=psw[c('eff0','eff1','effDiff')],
    SE=NA,
    psw=1
  )
}


proc1pstrata <- function(res1){
    if("pstrata"%in%names(res1))
        return(
            tibble(
                eff=c(0,1),
                estimator="pstrata",
                est=res1$pstrata[,"mean"],
                se=res1$pstrata[,"sd"],
                CIpercL=res1$pstrata[,"2.5%"],
                CIpercU=res1$pstrata[,"97.5%"]
            )
        ) else return(NULL)
}

proc1mest <- function(res1){
  tibble(
    eff = c(0,1),
    estimator= 'mest',
    est = res1$mest[c("eff0","eff1"),"estimates"],
    se = res1$mest[c("eff0","eff1"),"SE"],
  )
}

proc1pstrata <- function(res1){
  tibble(
    eff=c(0,1),
    estimator="pstrata",
    est=res1$pstrata[,"mean"],
    se=res1$pstrata[,"sd"],
    CIpercL=res1$pstrata[,"2.5%"],
    CIpercU=res1$pstrata[,"97.5%"]
  )
}

proc1bayes <- function(res1){
  tibble(
    eff = c(0,1),
    estimator= 'bayes',
    est = res1$bayes[c("eff0","eff1"),"mean"],
    se = res1$bayes[c("eff0","eff1"),"sd"],
    CIpercL = res1$bayes[c("eff0","eff1"), '2.5%'],
    CIpercU = res1$bayes[c("eff0","eff1"), '97.5%']
  )
}

proc1psw <- function(res1){
  if(is.null(res1$psw)) return(NULL)
  tibble(
    eff = c(0,1),
    estimator= 'psw',
    est = res1$psw[c("eff0","eff1")],
    se = NA
  )
}

proc1 <- function(res1){
  if(is.null(res1)) return(NULL)
  if(inherits(res1,'try-error')) return(NULL)

  res <- bind_rows(
    proc1mest(res1),
    proc1bayes(res1),
    proc1psw(res1),
    proc1pstrata(res1)
  )
  res$pop = vapply(res$eff, muEff, facs=res1$facs, numeric(1))
  res$samp <- unname(res1$true[paste0("S",res$eff)])
  res$rhat <- unname(res1$bayes[paste0("eff",res$eff),"Rhat"])
  res$auc <- attr(res1$mest, 'auc')

  res
}

proc <- function(res){
  out <- cbind(
    res[[1]]$facs,
    map_dfr(res,proc1)%>%
      mutate(
        CInormL = est - 2 * se,
        CInormU = est + 2 * se
      )
  )
  if("CIpercL"%in%names(out))
      out <- mutate(out,
        CIpercL = ifelse(is.na(CIpercL), CInormL, CIpercL),
        CIpercU = ifelse(is.na(CIpercU), CInormU, CIpercU)
      )
  out
}



loadRes <- function(ext1='',ext2='',pswResults){
    df <- FALSE
    if(file.exists(paste0('simResults',ext1,'/cases',ext2,'.RData'))){
        load(paste0('simResults',ext1,'/cases',ext2,'.RData'))
        df <- TRUE
    } else cases <- cbind(
               grep(paste0('sim[0-9]+',ext2,'\\.RData'),list.files(paste0("simResults",ext1,"/")),value=TRUE)
           )

    cat(nrow(cases),' cases to process\n')
    map_dfr(seq_len(nrow(cases)), function(i){
      cat(i,' ',sep='')
      load(if(df) paste0('simResults',ext1,'/sim',i,ext2,'.RData') else paste0("simResults",ext1,"/",cases[i,1]))
      resT <- proc(res)
      resT$run <- i
      resT})
}

bp <- function(pd,subset,facet,title=deparse(substitute(subset)),
                ylim,Labeller="label_value",labSize=5){


  r <- if (missing(subset))
    rep_len(TRUE, nrow(pd))
  else {
    e <- substitute(subset)
    r <- eval(e, pd, parent.frame())
    if (!is.logical(r))
      stop("'subset' must be logical")
    r & !is.na(r)
  }

  pd <- pd[r, , drop = TRUE]

  if(missing(ylim)) ylim = quantile(pd$errP, c(0.01, 0.99))

  pdSmall <- pd%>%
      group_by(across(all_of(rownames(attr(terms(facet),"factors")))))%>%
      sample_frac(0.1)%>%
      ungroup()

  p <- ggplot(pd,
         aes(
           x = estimator,
           y = errP,
           fill = estimator,
           color = estimator
         )) +
    geom_jitter(data=pdSmall,alpha = 0.1) +
    geom_violin(#position = "dodge2",
      color = 'black',
      quantile.linetype="solid",
      quantiles = 0.5) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim =ylim)+
    facet_nested(facet , #scales = "free",
               labeller = Labeller)+
    theme(legend.position = 'none')

    if(!is.null(title)) p <- p+ggtitle(title)

  if(!is.null(ylim))
    p=p+stat_summary(aes(estimator,errP),geom='label', fun.data=function(xx) data.frame(y=if(sum(xx<ylim[1])>0) ylim[1] else ylim[1]-100,label=paste('+',sum(xx <ylim[1]))),inherit.aes=FALSE,size=labSize)+
      stat_summary(aes(estimator,errP),geom='label', fun.data=function(xx) data.frame(y=if(sum(xx>ylim[2])>0) ylim[2] else ylim[2]+100,label=paste('+',sum(xx >ylim[2]))),inherit.aes=FALSE,size=labSize)
  p
}


load1 <- function(i,ext1='',ext2=''){
  load(paste0('simResults',ext1,'/sim',i,ext2,'.RData'))#fn[i]))
  list(res=res,facs=facs)
}

justLoad <- function(ext1='',ext2=''){
  load(paste0('simResults',ext1,'/cases',ext2,'.RData'))

  resList <- lapply(1:nrow(cases),load1)

  list(resList=resList,cases=cases)
}

justProc <- function(resList)
  map(1:length(resList),
      function(i){
        cat(i,'/',length(resList),' ')
        resT <- try(map_dfr(resList[[i]]$res,bayesProc,facs=resList[[i]]$facs))
        if(inherits(resT,'try-error')) resT <- as_tibble(resList[[i]]$facs)
        resT$run <- i
        resT
      })
