source("code/regression.r")
library(rstan)

psMod0 <- glm(S2~experiment_id+
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

summary(psMod0)
aucMod(psMod0)

psModStep <- step(psMod0)
summary(psModStep)
aucMod(psModStep)


effsFromFit(
  ps0 <- est(
    data = dat,
    covFormU = formula(psMod0)[-2],
    trt = "condition",
    out = "normalized_student_learning",
    use = "S2")
)

arm::binnedplot(predict(ps0$psMod,type="response"),
                resid(ps0$psMod,type="response"))
plot(ps0$outMod,which = 1)

fwlPlotDat <- function(mod,colNum){
  mf <- model.frame(mod)
  yName <- names(mf)[1]
  xName <- names(mf)[colNum]
  mY <- update(mod,as.formula(paste(".~.-",xName)),data=mf)
  mX <- update(mod,as.formula(paste(xName,"~.-",yName)),data=mf)
  sMod <- lm(resid(mY)~resid(mX))
  data.frame(resid=resid(sMod),fitted=fitted(sMod),regressor=xName)
  #ggplot(mapping=aes(x = fitted(sMod),
}

fwlPlots <- function(mod){
  mf <- model.frame(mod)
  which(vapply(mf[,-1],n_distinct,1)>10 &
          vapply(mf[,-1],is.numeric,TRUE))%>%
  map_dfr(~fwlPlotDat(mod,.+1))%>%
    ggplot(aes(fitted,resid))+geom_point()+geom_smooth()+
    facet_wrap(~regressor,scales="free_x")
}


fwlPlots(ps0$outMod)

effsFromFit(
  ps1 <- est(
    data = dat,
    covFormU = formula(psModStep)[-2],
    trt = "condition",
    out = "normalized_student_learning",
    use = "S2")
)

arm::binnedplot(predict(ps1$psMod,type="response"),
                resid(ps1$psMod,type="response"))

plot(ps1$outMod,which = 1)

effsFromFit(
  ps3 <- est(
    data = dat,
    covFormU = update(formula(psMod0)[-2],.~.
                      -student_prior_average_correctness
                      -class_prior_average_correctness),
    trt = "condition",
    out = "normalized_student_learning",
    use = "S2")
)


### psw
PSW <- dat%>%
  rename(Z=condition,
         Y=normalized_student_learning)%>%
  psw(psMod = psMod0)

### just to be sure, to GEEPERs in the same way
dat%>%
  rename(Z=condition,
         Y=normalized_student_learning)%>%
  est(psMod = psMod0)


### bayesian mixture model
Xout <- model.matrix(ps0$outMod)
Xout <- scale(Xout[,-c(1,which(colnames(Xout)%in%c('condition','Sp','condition:Sp')))])
Xu <- scale(model.matrix(formula(ps0$psMod)[-2],data=dat)[,-1])
sdat <- list(
  YctlY=dat$normalized_student_learning[dat$condition==0],
  YtrtY=dat$normalized_student_learning[dat$condition==1],
  XctlU=Xu[dat$condition==0,],
  XtrtU=Xu[dat$condition==1,],
  XtrtU=model.matrix(ps0$psMod)[,-1],
  XctlY=Xout[dat$condition==0,],
  XtrtY=Xout[dat$condition==1,],
  bottomOuter=dat$S2[dat$condition==1],
  nc=sum(1-dat$condition),
  nt=sum(dat$condition)
)
sdat$ncovU=ncol(sdat$XctlU)
sdat$ncovY=ncol(Xout)


psStan=stan('code/fh2t/psMod.stan',data=sdat,iter = 4000,cores = 4)
save(psStan,file='results/psEthanStan2.RData')

print(psStan,par=c("bottomOuterATE","notbottomOuterATE","ATEdiff"))


source("code/em.r")

print(em1 <- PS_EM(dat,psMod = psMod0,trt="condition",out="normalized_student_learning"))



### IV
library(AER)

iv <- ivreg(normalized_student_learning ~ S2 | condition, data = dat)
coeftest(iv,vcovHC)

iv2 <- ivreg(
    normalized_student_learning ~ S2 | condition + experiment_id +
    student_prior_started_skill_builder_count + student_prior_skill_builder_percent_completed +
    student_prior_started_problem_set_count + student_prior_problem_set_percent_completed +
    student_prior_completed_problem_count + student_prior_median_first_response_time +
    student_prior_median_time_on_task + student_prior_average_attempt_count +
    student_prior_average_correctness + class_age_in_days + class_student_count +
    class_prior_started_skill_builder_count + class_prior_skill_builder_percent_completed +
    class_prior_started_problem_set_count + class_prior_problem_set_percent_completed +
    class_prior_completed_problem_count + class_prior_median_first_response_time +
    class_prior_median_time_on_task + class_prior_average_attempt_count +
    class_prior_average_correctness + teacher_account_age_in_days +
    opportunity_zone,
    data = dat)

coeftest(iv2,vcovHC)
