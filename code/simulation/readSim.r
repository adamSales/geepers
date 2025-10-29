library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

source('code/simulation/readSimFuncs.r')

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

if(file.exists("simResults/lognormalResults.RData")&!fromScratch){
    load("simResults/lognormalResults.RData")
} else{
    resultsLognorm=loadRes(ext2='lognorm')
    save(resultsLognorm,file='simResults/lognormalResults.RData')
}


#resultsN100=loadRes(ext2='n100_mu01is0')
#save(resultsN100,file='simResults/resultsN100.RData')


#resultsNs=loadRes(ext2="ns")

#resultsB1s=loadRes(ext2="b1s")#s_mu01is0')

#results=loadRes()

#resB101=loadRes(ext2="b01_mu01is0")

### merge results by n
load("simResults/casesns.RData")#_mu01is0.RData')
casesNs=cases
#load('simResults/casesn100.RData')
#casesN100=cases
load('simResults/cases.RData')

resultsNs=bind_rows(
  resultsNs,
  # filter(resultsN100,mu01%in%resultsNs$mu01[1],
  #        errDist%in%resultsNs$errDist,
  #        b1%in%resultsNs$b1,
  #        intS%in%resultsNs$intS,
  #        intZ%in%resultsNs$intZ
  #        ),
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
