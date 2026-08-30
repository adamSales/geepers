
source('code/readSimFuncs.r')

if(!exists("fromScratch")) fromScratch <- TRUE
#### read, process results from main simulation
#print(load('simResults/pswResults.RData'))

if(file.exists("simResults/mainResults.RData")&!fromScratch){
    load("simResults/mainResults.RData")
} else{
    results=loadRes()
    save(results,file='simResults/mainResults.RData')
}

if(file.exists("simResults/nsResults.RData")&!fromScratch){
    load("simResults/nsResults.RData")
} else{
    resultsNs=loadRes(ext2='ns')
    save(resultsNs,file='simResults/nsResults.RData')
}

if(file.exists("simResults/b1sResults.RData")&!fromScratch){
    load("simResults/b1sResults.RData")
} else{
    resultsB1s=loadRes(ext2='b1s')
    save(resultsB1s,file='simResults/b1sResults.RData')
}


if(file.exists("simResults/lognormResults.RData")&!fromScratch){
    load("simResults/lognormResults.RData")
} else{
    resultsLN <- loadRes(ext2="lognorm")
    save(resultsLN,file="simResults/lognormResults.RData")
}

if(file.exists("simResults/oneVarIntResults.RData")&!fromScratch){
    load("simResults/oneVarIntResults.RData")
} else{
    resultsInt <- loadRes(ext2="oneVarInt")
    save(resultsInt,file="simResults/oneVarIntResults.RData")
}

if(file.exists("simResults/nox1Results.RData")&!fromScratch){
    load("simResults/nox1Results.RData")
} else{
    resultsNox1 <- loadRes(ext2="oneVarIntnox1")
    save(resultsNox1,file="simResults/nox1Results.RData")
}





results <- bind_rows(
    filter(results,errDist!="lognorm"),
    resultsLN)%>%
    filter(!intZ)%>%
    bind_rows(resultsInt)

save(results,file='simResults/mainResults.RData')


#resultsN100=loadRes(ext2='n100_mu01is0')
#save(resultsN100,file='simResults/resultsN100.RData')


### merge results by n
load("simResults/casesns.RData")#_mu01is0.RData')
casesNs=cases
#load('simResults/casesn100.RData')
#casesN100=cases
load('simResults/cases.RData')

resultsNs=bind_rows(
  resultsNs,

  filter(results,mu01%in%resultsNs$mu01[1],
         errDist%in%resultsNs$errDist,
         b1%in%resultsNs$b1,
         intS%in%resultsNs$intS,
         intZ%in%resultsNs$intZ
         ))%>%
    filter(estimator!="psw")

save(resultsNs,file='simResults/resultsNs.RData')

### merge results by B1s
casesTot=cases
load('simResults/casesb1s.RData')
casesB1=cases

resultsB1s=bind_rows(
    resultsB1s,
 #   resB101,
  filter(results,run%in%(
    left_join(casesB1[1,],casesTot%>%mutate(run=1:n())%>%select(-b1))%>%pull(run)
    ))
  )
save(resultsB1s,file='simResults/resultsB1s.RData')



### quick look
results%>%filter(b1>0,n==500,eff==1,errDist!="lognorm")%>%
    group_by(mu01,mu10,mu11,b1,errDist,intS,intZ,estimator)%>%
    summarize(
        bias=mean(est-pop),
        realSE=sd(est-samp),
        covr=mean(CIpercL<samp&CIpercU>samp))%>%print(n=Inf)


resultsLN%>%filter(b1>0,n==500,eff==1)%>%
    group_by(mu01,mu10,mu11,b1,errDist,intS,intZ,estimator)%>%
    summarize(
        bias=mean(est-pop),
        realSE=sd(est-samp),
        covr=mean(CIpercL<samp&CIpercU>samp))%>%print(n=Inf)


resultsLN%>%filter(b1>0,n==500,eff==1,estimator=="pstrata")%>%
    group_by(mu01,mu10,mu11,b1,errDist,intS,intZ,estimator)%>%
    summarize(
        bias=mean(est-pop),
        realSE=sd(est-samp),
        covr=mean(CIpercL<samp&CIpercU>samp))%>%print(n=Inf)
resultsLN%>%filter(b1>0,n==500,eff==1,estimator=="mest")%>%
    group_by(mu01,mu10,mu11,b1,errDist,intS,intZ,estimator)%>%
    summarize(
        bias=mean(est-pop),
        realSE=sd(est-samp),
        covr=mean(CIpercL<samp&CIpercU>samp))%>%print(n=Inf)


summ <- results%>%filter(b1>0,n==500,eff==1,intZ=TRUE)%>%
    group_by(mu01,mu10,mu11,b1,errDist,intS,intZ,estimator)%>%
    summarize(
        bias=mean(est-pop),
        realSE=sd(est-samp),
        covr=mean(CIpercL<samp&CIpercU>samp))%>%print(n=Inf)

filter(summ,intZ)%>%print(n=Inf)


resultsNox1%>%filter(estimator=="mest",b1>0,eff==1)%>%
    group_by(mu01,mu10,mu11,b1,errDist,intS,intZ,estimator)%>%
    summarize(
        bias=mean(est-pop),
        realSE=sd(est-samp),
        covr=mean(CInormL<samp&CInormU>samp))%>%print(n=Inf)
