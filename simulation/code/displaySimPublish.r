library(tidyverse)
library(ggh4x)
library(broom)
library(kableExtra)
library(gridExtra)
library(ggpubr)
library(tikzDevice)
select <- dplyr::select

options(
tikzLatexPackages = c(
    "\\usepackage{tikz}",
    "\\usepackage{bm}",
    "\\usepackage{amsmath}",
"\\usepackage[active,tightpage]{preview}",
"\\PreviewEnvironment{pgfpicture}",
"\\setlength\\PreviewBorder{0pt}"
),
tikzXelatexPackages = c(
"\\usepackage{tikz}\n",
"\\usepackage[active,tightpage,xetex]{preview}\n",
"\\usepackage{fontspec,xunicode}\n",
"\\PreviewEnvironment{pgfpicture}\n",
"\\setlength\\PreviewBorder{0pt}\n"
),
tikzLualatexPackages = c(
    "\\usepackage{tikz}\n",
    "\\usepackage{bm}",
        "\\usepackage{amsmath}",
"\\usepackage[active,tightpage,psfixbb]{preview}\n",
"\\usepackage{fontspec,xunicode}\n",
"\\PreviewEnvironment{pgfpicture}\n",
"\\setlength\\PreviewBorder{0pt}\n"
)
)

source('simulation/code/readSimFuncs.r')

#### after simulation results have been loaded and pre-processed...
## source('code/simulation/readSim.r')
load('simResults/mainResults.RData')
load('simResults/resultsNs.RData')
load('simResults/resultsB1s.RData')


###############################################################################
### results across ns
###############################################################################
 pdns <- resultsNs %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=c(0,0.2,0.5),
                labels=c(bquote(alpha==0),bquote(alpha==0.2),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction"),
    N=paste0("n=",n)
  )


### numbers of Rhat>1.1
nonconvN <- pdns%>%
    filter(estimator=="Mixture")%>%
    group_by(n,eff,estimator)%>%
    summarize(nonconv=sum(rhat>=1.1),percNonconv=nonconv/n())


#### figure in the paper
biasN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
    summarize(bias=mean(est-pop,na.rm=TRUE),
            ttest=tidy(t.test(est-pop)))%>%
  ungroup()%>%
  bind_cols(.$ttest)%>%
  select(-ttest)


seN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
  summarize(se=sd(errP,na.rm=TRUE),
            conf.high=sqrt((n()-1)/qchisq(0.025,n()-1))*se,
            conf.low=sqrt((n()-1)/qchisq(0.975,n()-1))*se)%>%ungroup()

#full_join(seN,biasN)%>%filter(eff==1,estimator!="PSW")%>%summarize(across(c(se,bias),max),coef=round(se/bias))

plotByN=bind_rows(
  seN%>%
  ## divide by 3 to plot w/bias; axis multiplied by 3 below
  mutate(across(c(se,conf.high,conf.low),~./3),meas="SE")%>%
  rename(what=se),
  biasN%>%mutate(meas="Bias")%>%rename(what=bias))%>%
    filter(eff==0,estimator!="PSW")%>%
    ggplot(aes(n,what,color=estimator,group=paste0(estimator,meas),linetype=meas,fill=estimator,shape=meas))+
    geom_errorbar(aes(ymin=conf.low,ymax=conf.high), width=10)+
    #geom_point()+
    geom_line()+geom_hline(yintercept=0)+#,linetype="dotted",size=2)+
    scale_y_continuous(name="Bias",sec.axis=sec_axis(trans=~.*3,name='Standard Error'))+
    scale_shape_manual(values=c(0,16))+
    annotate('text',600,.11,label=list(bquote(paste(#atop(
                                "Normal Resid., No Interactions, ",alpha==0.5))),#", "~mu[C]^1-mu[C]^0==0.3))),
             parse=TRUE)+
    scale_x_continuous(name="Sample Size Per Group (n)",breaks=c(seq(100,700,200),1000))+
  #ggtitle("Bias and Standard Error by n")+
    theme(legend.title=element_blank())

ggsave("paper/simFigs/biasSEbyN.jpg",plot=plotByN,width=6,height=3)




###############################################################################
### results across B1
###############################################################################

 pdb1 <- resultsB1s %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=seq(0,1,.1),
                labels=lapply(seq(0,1,.1), function(x) bquote(alpha==.(x)))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction"),
    N=paste0("n=",n)
  )

### numbers of Rhat>1.1
nonconvB1 <- pdb1%>%
    filter(estimator=="Mixture")%>%
    group_by(b1,eff,estimator)%>%
    summarize(nonconv=sum(rhat>=1.1),percNonconv=nonconv/n())



#### figure in the paper

biasB1 <- pdb1%>%
    filter(rhat<1.1)%>%
    group_by(b1,eff,estimator)%>%
    summarize(bias=mean(errP,na.rm=TRUE),
              ttest=tidy(t.test(errP)))%>%
    ungroup()%>%
    bind_cols(.$ttest)%>%
    select(-ttest)

  #summarize(bias=mean(errP,na.rm=TRUE))%>%ungroup()

seB1 <- pdb1%>%
    filter(rhat<1.1,!is.na(errP))%>%
    group_by(b1,eff,estimator)%>%
    summarize(se=sd(errP),
              bias=mean(errP),
              conf.high=sqrt((4999)/qchisq(0.025,4999))*se,
              conf.low=sqrt((4999)/qchisq(0.975,4999))*se)%>%ungroup()

#full_join(seB1,biasB1)%>%filter(eff==1,estimator!="PSW",b1>0)%>%summarize(across(c(se,bias),max),coef=round(se/bias))

plotByAlpha=bind_rows(
    seB1%>%mutate(across(c(se,conf.high,conf.low),~./5),meas="SE")%>%rename(what=se),
    biasB1%>%mutate(meas="Bias")%>%rename(what=bias))%>%
    filter(b1>0.1,eff==0,estimator!="PSW")%>%
    ggplot(aes(b1,what,color=estimator,group=paste0(estimator,meas),linetype=meas,fill=estimator,shape=meas))+
    #geom_point()+
    geom_line()+geom_hline(yintercept=0)+#,linetype="dashed")+
    geom_errorbar(aes(ymin=conf.low,ymax=conf.high), width=.01)+
    scale_y_continuous(name="Bias",sec.axis=sec_axis(trans=~.*5,name='Standard Error'))+
    scale_x_continuous(breaks=seq(0,1,0.2),minor_breaks=seq(0.1,0.9,0.2),labels=sprintf("%0.1f",seq(0,1,0.2)))+
    scale_shape_manual(values=c(0,16))+
    annotate('text',.7,.08,label=list(bquote(
                               paste(#atop(
                                   "Normal Resid., No Interactions, ",n==500))),# "~mu[C]^1-mu[C]^0==0.3))),
             parse=TRUE)+
    xlab(bquote(alpha))#+ggtitle(bquote("Bias and Standard Error by "~alpha))

ggsave("paper/simFigs/biasSEbyB1.jpg",plot=plotByAlpha,width=6,height=3)


#### combine results by n and B1
ggarrange(plotByN,plotByAlpha,ncol=1,common.legend = TRUE, legend = "bottom")
ggsave("paper/simFigs/biasSEbyB1n.jpg",width=5,height=4)






######################################################################
### main results

pd <- results %>%
    filter(estimator%in%c("bayes","mest","psw"))%>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - samp,
    B1 = factor(b1,levels=c(0,0.3,0.5),
                labels=paste0("$\\alpha=",c(0,0.3,0.5),"$")),
#                    c(bquote(alpha==0),bquote(alpha==0.3),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              ##labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
              labels=c(bquote(mu[c]^1==0),bquote(mu[c]^1==0.3))),
    PE=paste('Stratum',eff),
    dist=factor(paste(c(unif='Unif.',norm='Normal',lognorm="Lognorm.")[errDist],''),
                levels=paste(c("Normal","Lognorm.","Unif."),'')),
    #estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    estimator=c(bayes='\\textsc{pmm}',mest='\\textsc{geepers}',psw='\\textsc{psw}')[estimator],
    interactionZ=ifelse(intZ,"Z\ninteraction","No Z\ninteraction"),
    interactionS=ifelse(intS,"S\ninteraction","No S\ninteraction"),
    intAll=ifelse(intZ,ifelse(intS,"$\\bm{x}\\text{:}S_T$\\&\n$\\bm{x}\\text{:}Z$","$\\bm{x}\\text{:}Z$"),ifelse(intS,"$\\bm{x}\\text{:}S_T$","No\ninter.")),
    intAll=factor(intAll,levels=c("No\ninter.", "$\\bm{x}\\text{:}S_T$","$\\bm{x}\\text{:}Z$","$\\bm{x}\\text{:}S_T$\\&\n$\\bm{x}\\text{:}Z$"))
  ) #%>%
  #filter(estimator=='PSW'|rhat<1.1)


### numbers of Rhat>1.1
nonconvMain <- pd%>%
    filter(eff==1,mu01==0,n==500,estimator=='\\textsc{pmm}')%>%
    group_by(B1,dist,intAll,estimator,run)%>%
    summarize(nonconv=sum(rhat>=1.1),percNonconv=nonconv/n())


### figure for paper

violin <- bp(filter(pd,rhat<1.1),
   subset=b1 > 0&  eff==0&mu01==0&n==500&PE=='Stratum 0',
   title=NULL,#"Normal Residuals",
   ylim=c(-1.5,1.5),#c(-1.5,1.5),
   facet=B1~intAll+dist,#interactionZ+interactionS,
   Labeller=labeller(B1="none",#label_parsed,
                     intAll=label_value),labSize=2.5
   )+
    labs(y=bquote("Estimation Error for "~tau^0),#subtitle="No Interactions",
         x=NULL)+
  #theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
          #strip.text.y=element_blank(),strip.background.y=element_blank())+
                                        #plot.subtitle = element_text(size=10),)+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))





#pdf(
tikz("paper/simFigs/boxplotsNew.tex",#pdf",
     width=8,height=5,standAlone=TRUE)
print(violin)
                                        #grid.arrange(norm,unif,lognorm,nrow=1)
dev.off()

setwd("paper/simFigs")
system("lualatex boxplotsNew.tex")
setwd("../..")

#######################################################
#### rmse appendix
#######################################################

rmse <- pd%>%
  filter(rhat<1.1|estimator=='PSW')%>%
  group_by(n,mu01,errDist,b1,intS,intZ,interactionS,interactionZ,eff,estimator)%>%
   summarize(rmse=sqrt(mean(errP^2,na.rm=TRUE)))%>%ungroup()



sink('paper/rmseTabAppendix500.tex')

cbind(
rmse%>%filter(errDist!='mix',n==500,eff!="Diff",b1==.3,mu01==0)%>%
  transmute(`Residual\nDist.`=c(norm='Normal',lognorm="Lognormal",unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers",PSW="\\psw")[estimator],
#         `$\\beta_1$`=mu01,
         `Prin.\nEff`=paste0('$\\tau^',eff,'$'),
         rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse),
rmse%>%filter(errDist!='mix',n==500,mu01==0,eff!="Diff",b1==.5)%>%
  transmute(`Res.\nDist.`=c(norm='Norm.',unif='Unif.')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers",PSW="\\psw")[estimator],
         n,
         `$\\beta_1$`=mu01,
         `Prin.\nEff`=paste0('$\\tau^',eff,'$'),
         rmse,
         xx=rep(1:(n()/3),each=3))%>%
select(xx,estimator,rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse)%>%
select(-xx))%>%#GEEPERs,Mix.,PSW))%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 5,  "$\\\\alpha=0.3$" = 3, "$\\\\alpha=0.5$" = 3),escape=FALSE)%>%
     add_header_above(c(" " = 5,  "$n=500$"=6),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")

sink()




#######################################################
#### coverage
#######################################################
 coverage <- pd%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,interactionS,interactionZ,intS,intZ,eff,estimator)%>%
   summarize(coverage=mean(CInormL<=samp & CInormU>=samp,na.rm=TRUE))%>%ungroup()



redCov <- function(coverage) paste0("\\rd{",sprintf("%.2f",coverage),"}")
condRed <- function(coverage,intZ,intS,estimator,errDist,b1=1)
    ifelse(intZ|intS,redCov(coverage),
    ifelse(errDist%in%c("unif","lognorm")&estimator=="\\textsc{pmm}",redCov(coverage),
    ifelse(estimator=="\\textsc{geepers}"&b1==0,redCov(coverage),sprintf("%.2f",coverage))))


sink('paper/coverageTab.tex')
map(c(0,0.3,0.5),function(b){
   tab <- coverage%>%
      mutate(coverage=condRed(coverage,intZ,intS,estimator,errDist,b1))%>%
      filter(n==500,eff==0,errDist!='mix',mu01==0,b1==b,estimator!="\\textsc{psw}")%>%
       transmute(`Residual\nDist.`=
                      factor(c(norm='Normal',unif='Uniform',lognorm="Lognormal")[errDist],
                            levels=c("Normal","Lognormal","Uniform")),
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers")[estimator],
         coverage)%>%
       pivot_wider(names_from=estimator,values_from=coverage)%>%
       arrange(`Residual\nDist.`)
   if(b==0) return(tab) else return(tab[,4:5])
}
)%>%
    bind_cols(.name_repair="minimal")%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 3, "$\\\\alpha=0$" = 2,
                       "$\\\\alpha=0.3$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    collapse_rows(columns=1,latex_hline="major",valign="middle")%>%
    footnote(general=c("\\\\footnotesize Based on 5000 replications. $n=500$. Simulation standard error $\\\\approx 1/3$ percentage point. Estimates colored \\\\rd{red} indicate cases where the assumptions of the model are not met."),escape=FALSE,footnote_as_chunk = TRUE,threeparttable=TRUE)%>%print()
sink()


#######################################################
#### type-1 error
#######################################################
type1err <- pd%>%
  filter(rhat<1.1,estimator!='\\textsc{psw}',eff==0,!intZ)%>%
  group_by(n,mu01,errDist,dist,b1,intS,intZ,interactionS,interactionZ,eff,estimator)%>%
   summarize(type1err=mean(abs(est)>2*se))%>%ungroup()


sink("paper/typeIerror.tex")
map(c(0,0.3,0.5),function(b){
    tab <- type1err%>%
      mutate(type1err=condRed(type1err,intZ,intS,estimator,errDist,b1))%>%
      filter(n==500,eff==0,errDist!='mix',mu01==0,b1==b,estimator!="\\textsc{psw}")%>%
       transmute(`Residual\nDist.`=
                      factor(c(norm='Normal',unif='Uniform',lognorm="Lognormal")[errDist],
                            levels=c("Normal","Lognormal","Uniform")),
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers")[estimator],
         type1err)%>%
       pivot_wider(names_from=estimator,values_from=type1err)%>%
       arrange(`Residual\nDist.`)
   if(b==0) return(tab) else return(tab[,4:5])
})%>%
    bind_cols(.name_repair="minimal")%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 3, "$\\\\alpha=0$" = 2,
                       "$\\\\alpha=0.3$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    collapse_rows(columns=1,latex_hline="major",valign="middle")%>%
    footnote(general=c("\\\\footnotesize Based on 5000 replications. $n=500$. Simulation standard error $\\\\approx 1/3$ percentage point. Estimates colored \\\\rd{red} indicate cases where the assumptions of the model are not met."),escape=FALSE,footnote_as_chunk = TRUE,threeparttable=TRUE)%>%print()
sink()

########################
## coverage tab for appendix: n=500
########################

sink('paper/coverageTabAppendix500.tex')
cbind(
    coverage%>%filter(errDist!='mix',eff!="Diff",b1==0,estimator!='\\textsc{psw}',n==500
                      )%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         `$\\beta_1$`=mu01,
         `Prin.\nEff`=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage),
coverage%>%filter(errDist!='mix',eff!="Diff",b1==.3,estimator!='PSW',n==500)%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         `$\\beta_1$`=mu01,
         `Prin.\nEff`=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`\\textsc{geepers}`,`\\textsc{pmm}`),
coverage%>%filter(errDist!='mix',n==500,eff!="Diff",b1==.5,estimator!='\\textsc{psw}')%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         n,
         `$\\beta_1$`=mu01,
         `Prin.\nEff`=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`\\textsc{geepers}`,`\\textsc{pmm}`))%>%
kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 5, "$\\\\alpha=0$" = 2, "$\\\\alpha=0.3$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    add_header_above(c(" " = 5, "$n=500$"=6),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")
sink()


### how does AUC very w b1 and n?
resultsB1s%>%
  group_by(b1)%>%
  mutate(meanAUC=mean(auc))%>%
  ggplot(aes(as.factor(b1),auc))+geom_boxplot()+geom_point(aes(y=meanAUC))+geom_smooth(se=FALSE)+ylab('AUC')+xlab(bquote(alpha))
  ggsave('paper/simFigs/alphaAUC.pdf',width=5,height=4)
                                        #  scale_y_continuous('Avg. AUC',seq(.5,1,.1))




############################################################################
## no x1
############################################################################

load("simResults/nox1Results.RData")


pdnoX1 <- resultsNox1 %>%
    filter(estimator%in%c("bayes","mest","psw"))%>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - samp,
    B1 = factor(b1,levels=c(0,0.3,0.5),
                labels=paste0("$\\alpha=",c(0,0.3,0.5),"$")),
#                    c(bquote(alpha==0),bquote(alpha==0.3),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              ##labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
              labels=c(bquote(mu[c]^1==0),bquote(mu[c]^1==0.3))),
    PE=paste('Stratum',eff),
    dist=factor(paste(c(unif='Unform',norm='Normal',lognorm="Lognormal")[errDist],'Errors'),
                levels=paste(c("Normal","Lognormal","Unform"),'Errors')),
    #estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    estimator=c(bayes='\\textsc{pmm}',mest='\\textsc{geepers}',psw='\\textsc{psw}')[estimator],
    interactionZ=ifelse(intZ,"Z\ninteraction","No Z\ninteraction"),
    interactionS=ifelse(intS,"S\ninteraction","No S\ninteraction"),
    intAll=ifelse(intZ,ifelse(intS,"$\\bm{x}\\text{:}S_T$\\&\n$\\bm{x}\\text{:}Z$","$\\bm{x}\\text{:}Z$"),ifelse(intS,"$\\bm{x}\\text{:}S_T$","No\ninter.")),
    intAll=factor(intAll,levels=c("No\ninter.","$\\bm{x}\\text{:}Z$", "$\\bm{x}\\text{:}S_T$","$\\bm{x}\\text{:}S_T$\\&\n$\\bm{x}\\text{:}Z$"))
  ) #%>%
  #filter(estimator=='PSW'|rhat<1.1)


### figure for paper
violinNoX1 <- bp(pdnoX1,
   subset=b1 > 0&  eff==0&mu01==0&n==500&PE=='Stratum 0',
   title=NULL,#"Normal Residuals",
   ylim=c(-1.5,1.5),#c(-1.5,1.5),
   facet=B1~dist+intAll,#interactionZ+interactionS,
   Labeller=labeller(B1="none",#label_parsed,
                     intAll=label_value),labSize=2.5
   )+
    labs(y=bquote("Estimation Error for "~tau^1),#subtitle="No Interactions",
         x=NULL)+
  #theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
          #strip.text.y=element_blank(),strip.background.y=element_blank())+
                                        #plot.subtitle = element_text(size=10),)+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))


#pdf(
tikz("paper/simFigs/boxplotsNoX1.tex",#pdf",
     width=7,height=4,standAlone=TRUE)
print(violinNoX1)
                                        #grid.arrange(norm,unif,lognorm,nrow=1)
dev.off()

setwd("paper/simFigs")
system("lualatex boxplotsNoX1.tex")
setwd("../..")

coverageNoX1 <- pdnoX1%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,interactionS,interactionZ,intS,intZ,eff,estimator)%>%
   summarize(coverage=mean(CInormL<=samp & CInormU>=samp,na.rm=TRUE))%>%ungroup()




sink('paper/coverageTabNoX1.tex')

map(c(0.3,0.5),function(b){
   tab <- coverageNoX1%>%
      mutate(coverage=condRed(coverage,intZ,intS,estimator,errDist,b1))%>%
      filter(n==500,eff==0,errDist!='mix',mu01==0,b1==b,estimator!="\\textsc{psw}")%>%
       transmute(`Residual\nDist.`=
                      factor(c(norm='Normal',unif='Uniform',lognorm="Lognormal")[errDist],
                            levels=c("Normal","Lognormal","Uniform")),
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers")[estimator],
         coverage)%>%
       pivot_wider(names_from=estimator,values_from=coverage)
   if(b==0.3) return(tab) else return(tab[,4:5])
}
)%>%
    bind_cols(.name_repair="minimal")%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 3, "$\\\\alpha=0$" = 2,
                       "$\\\\alpha=0.3$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    collapse_rows(columns=1,latex_hline="major",valign="middle")%>%
    footnote(general=c("\\\\footnotesize Based on 5000 replications. $n=500$. Simulation standard error $\\\\approx 1/3$ percentage point. Estimates colored \\\\rd{red} indicate cases where the assumptions of the model are not met."),escape=FALSE,footnote_as_chunk = TRUE,threeparttable=TRUE)%>%print()
sink()




biasNoX1 <- pdnoX1%>%
    filter(eff==0,b1>0,!intS)%>%
    group_by(n,mu01,errDist,b1,interactionS,interactionZ,intS,intZ,eff,estimator)%>%
    summarize(bias=mean(est-samp),seBias=sd(est-samp)/sqrt(n()))



biasTab <-
    inner_join(
        biasNoX1%>%ungroup()%>%select(n:b1,intS:seBias),
        pd%>%
        filter(eff==0,b1>0,!intS)%>%
        group_by(n,mu01,errDist,b1,interactionS,interactionZ,intS,intZ,eff,estimator)%>%
        summarize(bias2=mean(est-samp),seBias2=sd(est-samp)/sqrt(n()))%>%ungroup()%>%select(n:b1,intS:seBias2))#%>%
#    mutate(bias=paste0(sprintf("%.2f", bias)," (",sprintf("%.2f", seBias),")"))

coverageTab <-
    inner_join(
        coverageNoX1%>%ungroup()%>%select(n:b1,intS:coverage)%>%filter(eff==0,b1>0,!intS)%>%rename(bias=coverage),
        coverage%>%ungroup()%>%select(n:b1,intS:coverage)%>%filter(eff==0,b1>0,!intS)%>%rename(bias2=coverage)
        )%>%
    filter(estimator!="\\textsc{psw}")%>%
    mutate(what="Coverage")%>%select(what,everything())




sink("paper/biasCoverageWithWithoutX1.tex")
bind_rows(
    inner_join(
        filter(biasTab,errDist=="norm")%>%select(b1,estimator,bias,bias2),
        filter(biasTab,errDist=="lognorm")%>%select(b1,estimator,bias,bias2),
        by=c("b1", "estimator"))%>%
    inner_join(
        filter(biasTab,errDist=="unif")%>%select(b1,estimator,bias,bias2),
        by=c("b1", "estimator"))%>%
    mutate(what="Bias")%>%select(what,everything()),
    inner_join(
        filter(coverageTab,errDist=="norm")%>%select(what,b1,estimator,bias,bias2),
        filter(coverageTab,errDist=="lognorm")%>%select(what,b1,estimator,bias,bias2),
        by=c("b1", "estimator","what"))%>%
    inner_join(
        filter(coverageTab,errDist=="unif")%>%select(what,b1,estimator,bias,bias2),
        by=c("b1", "estimator","what")))%>%
    rename("$\\alpha$"=b1)%>%
    kbl('latex',booktabs=TRUE,col.names=c(" ",names(.)[2:3],rep(c("Just $x_2$","$x_1$ \\& $x_2$"),3)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 3, "Normal"=2,"Lognormal"=2,"Uniform"=2))%>%
    add_header_above(c(" " = 3, "Error Distribution"=6))%>%
    collapse_rows(columns=1,latex_hline="major",valign="middle")%>%
    collapse_rows(columns=2,latex_hline="none",valign="middle")%>%
    footnote(general=c("\\\\footnotesize Based on 5000 replications. $n=500$. Simulation standard error $\\\\approx 1/3$ percentage point."),escape=FALSE,footnote_as_chunk = TRUE,threeparttable=TRUE)%>%print()
sink()

### standard error estimates?
results%>%filter(estimator=="mest",eff==1,b1>0,se<100)%>%group_by(run,n,b1,errDist,intS,intZ)%>%
  summarise(trueVar=var(est),EestVar=mean(se^2),seError=EestVar-trueVar,
            pval=t.test(se^2,mu=trueVar)$p.value)%>%arrange(seError)%>%print(n=Inf)


results%>%filter(estimator=="mest",eff==0,b1>0,se<100)%>%group_by(run,n,b1,errDist,intS,intZ)%>%
  summarise(trueVar=var(est),EestVar=mean(se^2),seError=EestVar-trueVar,
            pval=t.test(se^2,mu=trueVar)$p.value)%>%arrange(seError)%>%print(n=Inf)

results%>%filter(estimator=="mest",eff==1,b1>0,run==14)%>%pull(se)%>%hist()


results%>%filter(estimator=="bayes",eff==1,b1>0,se<100)%>%group_by(run,n,b1,errDist,intS,intZ)%>%
  summarise(trueVar=var(est),EestVar=mean(se^2),seError=EestVar-trueVar,
            pval=t.test(se^2,mu=trueVar)$p.value)%>%arrange(seError)%>%print(n=Inf)


results%>%filter(estimator=="mest",eff==0,b1>0,se<100)%>%group_by(run,n,b1,errDist,intS,intZ)%>%
  summarise(trueVar=var(est),EestVar=mean(se^2),seError=EestVar-trueVar,
            pval=t.test(se^2,mu=trueVar)$p.value)%>%arrange(seError)%>%print(n=Inf)


resultsNs%>%filter(estimator=="mest",eff==1,b1>0,se<100)%>%group_by(run,n,b1,errDist,intS,intZ)%>%
  summarise(trueVar=var(est),EestVar=mean(se^2),seError=EestVar-trueVar,
            pval=t.test(se^2,mu=trueVar)$p.value)%>%arrange(n)%>%print(n=Inf)

resultsB1s%>%filter(estimator=="mest",eff==1,b1>0,se<100)%>%group_by(run,n,b1,errDist,intS,intZ)%>%
  summarise(trueVar=var(est),EestVar=mean(se^2),seError=EestVar-trueVar,
            pval=t.test(se^2,mu=trueVar)$p.value)%>%arrange(b1)%>%print(n=Inf)
