
dat <- read_csv("data/as.csv")

print("covariates:")
print(covVariables <- c(grep("_prior_",names(dat),value=TRUE),
                  "class_age_in_days","class_student_count","opportunity_zone"))


trans <- function(x){
    x <- ifelse(is.na(x),mean(x,na.rm=TRUE),x)
    if(is.numeric(x) & n_distinct(x)>2) x <- (x-mean(x))/sd(x)
    x
}


dat <- dat%>%
  mutate(
    opportunity_zone=opportunity_zone=="Yes",
    across(all_of(covVariables),trans))




psMod0 <- glm(S~experiment_id+
                #no_student_prior_skill_builders+
                #no_student_prior_problem_sets+no_student_prior_attempted_problems+
                #no_student_prior_completed_problems+
                student_prior_started_skill_builder_count+
                student_prior_skill_builder_percent_completed+
                student_prior_started_problem_set_count+
                student_prior_problem_set_percent_completed+
                student_prior_completed_problem_count+
                student_prior_median_first_response_time+
                student_prior_median_time_on_task+
                student_prior_average_attempt_count+
                student_prior_average_correctness+
                class_age_in_days+
                class_student_count+
                #no_class_prior_skill_builders+no_class_prior_problem_sets+
                #no_class_prior_attempted_problems+no_class_prior_completed_problems+
                class_prior_started_skill_builder_count+class_prior_skill_builder_percent_completed+
                class_prior_started_problem_set_count+class_prior_problem_set_percent_completed+
                class_prior_completed_problem_count+class_prior_median_first_response_time+
                class_prior_median_time_on_task+class_prior_average_attempt_count+
                class_prior_average_correctness+teacher_account_age_in_days+
                opportunity_zone,
             data=dat,
             subset=condition==1,
             family=binomial)


ate <- with(t.test(normalized_student_learning~condition,data=dat),
            data.frame(estimates=unname(estimate[2]-estimate[1]),
                       se=stderr,
                       ymin=-conf.int[2],
                       ymax=-conf.int[1]))



psModStep <- step(psMod0,trace=0)



### compare models
psMods <- c(
    "psMod0",
    "psModStep"
)%>%
    setNames(.,.)%>%
    map(get)


sapply(psMods,predict,type="response")%>%bind_cols()%>%pairs()
sapply(psMods,predict,type="response")%>%bind_cols()%>%cor()


do.call("anova",unname(psMods))


sapply(psMods,AIC)

sapply(psMods,aucMod)

set.seed(248)
nt <- sum(dat$condition)
folds <- sample(rep(1:10,ceiling(nt))[1:nt])
cv <- function(mod){
    pred.out <- numeric(sum(dat$condition))
    for(ff in 1:10){
        modStar <- update(mod,data=subset(dat,condition==1)[folds!=ff,])
        pred.out[folds==ff] <- predict(modStar,subset(dat,condition==1)[folds==ff,])
    }
    pred.out
}

newPreds <- sapply(psMods,cv)

apply(newPreds,2,function(x) auc(x,1-psMod0$y))



### binned residual plots
par(mfrow=c(2,1))
map(names(psMods),~arm::binnedplot(fitted(psMods[[.]]),resid(psMods[[.]],type="response"),main=.))

map(names(psMods),~arm::binnedplot(plogis(newPreds[,.]),psMod0$y-plogis(newPreds[,.]),main=.))


###### initial estimates
print(
    ests1 <- map(
        psMods,#[-1],
        ~est(data=dat,
             covFormU=formula(.)[-2],
             trt = "condition",
             out = "normalized_student_learning",
             use = "S")
    )
)



#####
### removing certain variables

rmGrepForm <- function(formula,pattern,verbose=TRUE,...){
    ttt <- terms(formula)
    predictors <- attr(ttt,"term.labels")
    oneside <- attr(ttt,"response")==0

    rem <- grep(pattern,predictors,...)
    if(verbose) cat("removing: ",predictors[rem],"\n")
    predictors <- predictors[-rem]

    upForm <- paste(
        if(!oneside) "." else "",
        "~",
        paste(predictors,collapse="+"))

    update(formula,upForm)
}

### remove prior correctness
effsFromFit(
  ests1$noPretest <- est(
    data = dat,
    covFormU = rmGrepForm(formula(psModStep)[-2],"correctness"),
    trt = "condition",
    out = "normalized_student_learning",
    use = "S")
)

### remove experiment indicators
effsFromFit(
  ests1$noExp <- est(
    data = dat,
    covFormU = update(formula(psModStep)[-2],~.
                      -experiment_id),
    trt = "condition",
    out = "normalized_student_learning",
    use = "S")
)

### remove completenss
print(
    ests1$noComplete <-
    est(
    data = dat,
    covFormU = rmGrepForm(formula(psModStep)[-2],"compl"),
    trt = "condition",
    out = "normalized_student_learning",
    use = "S")
)

cvAUC <- function(mod){
    pred.out <- numeric(sum(dat$condition))
    form <- formula(mod)
    for(ff in 1:10){
        modStar <- glm(formula(mod),family=binomial,data=model.frame(mod)[folds!=ff,])
        pred.out[folds==ff] <- predict(modStar,model.frame(mod)[folds==ff,])
    }
    auc(pred.out,1-mod$y)
}

cvAUCs <- map_dbl(ests1,~cvAUC(.$psMod))


save(ests1,ate,cvAUCs,file="results/ests1.RData")

#######################################################
### other methods
#######################################################

### psw
pswResults <-
    dat%>%
    rename(Z=condition,
           Y=normalized_student_learning)%>%
    psw(psMod = psModStep)



### bayesian mixture model
Xu <- Xout <- #scale(
          model.matrix(formula(psModStep)[-2],data=dat)[,-1]#)

sdat <- list(
  YctlY=dat$normalized_student_learning[dat$condition==0],
  YtrtY=dat$normalized_student_learning[dat$condition==1],
  XctlU=Xu[dat$condition==0,],
  XtrtU=Xu[dat$condition==1,],
  XctlY=Xout[dat$condition==0,],
  XtrtY=Xout[dat$condition==1,],
  complier=dat$S[dat$condition==1],
  nc=sum(1-dat$condition),
  nt=sum(dat$condition)
)
sdat$ncovU=ncol(sdat$XctlU)
sdat$ncovY=ncol(Xout)


psStan=stan('code/psMod.stan',data=sdat,iter = 4000,cores = 4)
sims <- rstan::extract(psStan)
save(psStan,file='results/psEthanStan.RData')



pmmResults <- summary(psStan,par=c("complierATE","neverTakerATE","ATEdiff"))$summary

draws <- sample(1:length(sims[["alphaTBO"]]),100)

predYt1 <- function(i,X,S,alphaTBO,alphaTNBO,betaY,sigTBO,sigTNBO)
    rnorm(1,crossprod(X[i,],betaY)+ifelse(S[i]==1,alphaTBO,alphaTNBO),
          ifelse(S[i]==1,sigTBO,sigTBO))

ytpred <- matrix(nrow=sdat$nt,ncol=100)
for(k in 1:100)
    ytpred[,k] <- map_dbl(1:sdat$nt,
                          ~predYt1(i=.x,
                                   X=sdat$XtrtY,
                                   S=sdat$complier,
                                   alphaTBO=sims$alphaTBO[k],
                                   alphaTNBO=sims$alphaTNBO[k],
                                   betaY=sims$betaY[k,],
                                   sigTBO=sims$sigTBO[k],
                                   sigTNBO=sims$sigTNBO[k]))

predYc1 <- function(i,X,piC,alphaCBO,alphaCNBO,betaY,sigTBO,sigTNBO){
    S <- rbinom(1,1,piC[i])
    rnorm(1,crossprod(X[i,],betaY)+ifelse(S==1,alphaCBO,alphaCNBO),
          ifelse(S==1,sigTBO,sigTBO))
}

plot(density(sdat$YtrtY))
for(k in 1:100) lines(density(ytpred[,k]),col="red")
lines(density(sdat$YtrtY),lwd=2)

ycpred <- matrix(nrow=sdat$nc,ncol=100)
for(k in 1:100)
    ycpred[,k] <- map_dbl(1:sdat$nc,
                          ~predYc1(i=.x,
                                   X=sdat$XctlY,
                                   piC=sims$piC[k,],
                                   alphaCBO=sims$alphaCBO[k],
                                   alphaCNBO=sims$alphaCNBO[k],
                                   betaY=sims$betaY[k,],
                                   sigTBO=sims$sigTBO[k],
                                   sigTNBO=sims$sigTNBO[k]))


plot(density(sdat$YctlY))
for(k in 1:100) lines(density(ycpred[,k]),col="red")
lines(density(sdat$YctlY),lwd=2)


##### IV
iv <- ivreg(normalized_student_learning ~ S | condition, data = dat)
coeftest(iv,vcovHC)


ivForm <- as.formula(paste0("normalized_student_learning ~ S | condition+",
                            as.character(formula(psModStep))[3]))

iv2 <- ivreg(ivForm,data=dat)


ivResults <- coeftest(iv2,vcovHC)

save(pswResults,pmmResults,ivResults,iv2,file="results/competingMethods.RData")

