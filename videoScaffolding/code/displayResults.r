
dat <- read_csv("data/as.csv")

expDesc <- read_csv("rawdata/experiment_descriptions.csv")

expDesc <- expDesc[expDesc$`Experiment ID`%in%dat$experiment_id,]



expTab <- dat%>%
    bind_rows(mutate(dat,experiment_id="tot"))%>%
    group_by(experiment_id)%>%
    summarize(Trt=sum(condition),Ctl=n()-Trt,`$\\bar{S}|Z=1$`=(sum(S)/Trt),
              `Avg. \\# Scaffolds|S=1`=mean(scaffold_problems_given[S==1],na.rm=T),
              yt=mean(normalized_student_learning[condition==1]),
              yc=mean(normalized_student_learning[condition==0]))%>%
    #inner_join(select(expDesc,experiment_id=`Experiment ID`,`Dependent Variable`))%>%
    mutate(Experiment=c(LETTERS[1:n()-1],"Total"))%>%select(Experiment, everything())

## %>%
##     bind_rows(
##         dat%>%summarize(Trt=sum(condition),Ctl=n()-Trt,`$\\bar{S}|Z=1$`=sum(S),
##               `Avg. # Given|S=1`=mean(scaffold_problems_given[S==1],na.rm=T))%>%
##     #inner_join(select(expDesc,experiment_id=`Experiment ID`,`Dependent Variable`))%>%
##     mutate(Experiment="Total")%>%select(Experiment, everything())
##     )

sink("results/expTab.tex")
expTab%>%
    select(-experiment_id)%>%
     kbl('latex',booktabs=TRUE,escape=FALSE,digits=2,
         col.names=c(names(.)[1:5],"Trt","Ctl"))%>%
    add_header_above(c(" " = 1, "Sample Size" = 2,"Compliance"=2,"Mastery Speed"=2))
sink()


                                        #,-experiment_id)


dat%>%filter(condition==1)%>%
    group_by(experiment_id,scaffold_problems_given)%>%
    summarize(n=n())%>%
    ungroup()%>%
    mutate(Experiment=expTab$Experiment[match(experiment_id,expTab$experiment_id)])%>%
    ggplot(aes(scaffold_problems_given,n))+geom_col()+
    facet_wrap(~Experiment,scales="free_y",ncol=1)

load("results/ests1.RData")

### plot results
modelNames=c(
    psMod0="Full",
    psModStep="AIC Opt",
    noPretest="-Pretest",
    noExp="-Exp. IDs",
    noComplete="-Completion",
    psw="\\textsc{psw}",
    pmm="\\textsc{pmm}",
    iv="\\textsc{iv}",
    ATE="")

geepersEstimates <-
    map(ests1,effsFromFit)%>%
    map_dfr(~tibble(
                eff=rownames(.),
                estimates=.[,"estimates"],
                se=.[,"SE"],
                ymin=estimates-2*.[,"SE"],
                ymax=estimates+2*.[,"SE"]#,
#                AUC=attr(.,"auc")
                ))%>%
    bind_cols(data.frame(model=rep(names(ests1),each=3)))%>%
    mutate(AUC=cvAUCs[model])%>%
    bind_rows(mutate(ate,model="ATE",eff="ATE"))%>%
    filter(eff!="diff")

geepersEstimates%>%
    mutate(
        eff=factor(c(eff0="Never-Takers",eff1="Compliers",ATE="ATE")[eff],
                   levels=c("Never-Takers","ATE","Compliers")),
        model=factor(modelNames[model],levels=modelNames)
    )%>%
    ggplot(aes(model,estimates,ymin=ymin,ymax=ymax))+
    geom_point()+
    geom_errorbar(width=0,linewidth=1)+
    facet_wrap(~eff,nrow=1,scales="free_x",space="free_x")+
    geom_hline(yintercept=0)+geom_hline(yintercept=ate$estimates,linetype="dotted")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    xlab(NULL)+ylab("Principal Effect")
ggsave("figure/prinEffs.jpg",width=6,height=2.5,units="in")


## alternative methods
print(load("results/competingMethods.RData"))

altEstimates <- bind_rows(
     with(as.data.frame(pswResults$coef),
         tibble(
             model="psw",
             eff=rownames(pswResults$coef),
             estimates=est,
             se=se,
             ymin=estimates-2*se,
             ymax=estimates+2*se)),
    with(as.data.frame(pmmResults),
          tibble(
              model="pmm",
              eff=c(complierATE="eff1",
                    neverTakerATE="eff0",
                    ATEdiff="diff")[rownames(pmmResults)],
              estimates=mean,
              se=sd,
              ymin=`2.5%`,
              ymax=`97.5%`)),
    tibble(
        model="iv",
        eff=c("eff0","eff1"),
        estimates=c(0,coef(iv2)["STRUE"]),
        se=c(0,coeftest(iv2,vcov.=vcovHC)["STRUE","Std. Error"]),
        ymin=c(0,coefci(iv2,vcov.=vcovHC)["STRUE","2.5 %"]),
        ymax=c(0,coefci(iv2,vcov.=vcovHC)["STRUE","97.5 %"]))
)%>%
    filter(!grepl("diff",eff))

tikz("figure/compareMethods.tex",width=6,height=2.5,standAlone=TRUE)
print(
    bind_rows(
        filter(geepersEstimates,model%in%c("psModStep","ATE")),
        altEstimates)%>%
    mutate(
        eff=factor(c(eff0="Never-Takers",eff1="Compliers",ATE="ATE")[eff],
                   levels=c("Never-Takers","ATE","Compliers")),
        model=factor(
            ifelse(model=="psModStep","\\textsc{geepers}",modelNames[model]),
            levels=c("\\textsc{geepers}",modelNames)))%>%
    ggplot(aes(model,estimates,ymin=ymin,ymax=ymax))+
    geom_point()+
    geom_errorbar(width=0,linewidth=1)+
    facet_wrap(~eff,nrow=1,scales="free_x",space="free_x")+
    geom_hline(yintercept=0)+geom_hline(yintercept=ate$estimates,linetype="dotted")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    xlab(NULL)+ylab("Principal Effect")
)
dev.off()

setwd("figure")
system("lualatex compareMethods.tex")
setwd("..")



bind_rows(
    mutate(
        filter(geepersEstimates,model!="ATE"),Method="\\textsc{geepers}"),
    mutate(
        altEstimates,
        Method=paste0("\\textsc{",model,"}"),
        model="psModStep")
)%>%
    mutate(
        eff=factor(c(eff0="Never-Takers",eff1="Compliers",ATE="ATE")[eff],
                   levels=c("Never-Takers","ATE","Compliers")),
        `PS Model`=factor(modelNames[model],levels=modelNames),
        Estimate=paste0(
            sprintf("%.2f",round(estimates,2)),
            " (",sprintf("%.2f",round(se,2)),")"),
        AUC=ifelse(is.na(AUC),"",sprintf("%.2f",round(AUC,2)))
)%>%
    pivot_wider(
        id_cols=c(Method,`PS Model`,AUC),names_from=eff,values_from=Estimate)%>%
    kbl(format="latex",booktabs = T, align = "c",escape=FALSE) %>%
    collapse_rows(columns = 1:2, latex_hline = "major", row_group_label_position = "first")%>%
    add_header_above(c(" "=3,"Principal Effects"=2))%>%
    gsub("\\\\bottomrule",
         paste0(
             "\\\\cmidrule{1-5} \n \\\\textsc{ATE} & & & \\\\multicolumn{2}{c}{",
             sprintf("%.2f",round(ate$estimates,2))," (",
             sprintf("%.2f",round(ate$se,2)),")}\\\\\\\\\n\\\\bottomrule"),.)%>%
    cat(file="results/estimatesTab.tex")


ciFunc <- function(est,ymin){
    me <- est-ymin
    paste0(
        sprintf("%.2f",round(est,2)),
        "$\\pm$",
        sprintf("%.2f",round(me,2)))
}

geepersEstimates%>%
    arrange(eff)%>%
    mutate(ci=ciFunc(estimates,ymin))

altEstimates%>%
    arrange(eff)%>%
    mutate(ci=ciFunc(estimates,ymin))


### diagnostic plots for appendix
pdf("figure/binnedResids.pdf",width=6.5,height=7.5)

binnedResids <- function(mod,...)
    arm::binnedplot(fitted(mod),resid(mod,type="response"),...)

par(mfrow=c(3,2))
map(names(ests1),~binnedResids(ests1[[.]]$psMod,main=modelNames[.]))
dev.off()


pdf("figure/outResids.pdf",width=6.5,height=7.5)
par(mfrow=c(3,2))
map(names(ests1),~plot(ests1[[.]]$outMod,which=1,main=modelNames[.]))
dev.off()



### coefficients for appendix
abbv <- c(student="stud.",average="avg",median="med.",percent="pct.",completed="comp.",
          response="resp.",correctness="crct.")
tr <- ests1%>%map(~.$psMod)%>%texreg(single.row=TRUE,stars=numeric(0),caption="Coefficient estimates for the principal score logistic regressions described in Section \\ref{sec:application}",label="tab:psMods") ## put in texreg options
for(ex in expTab$experiment_id[1:3]) tr <- gsub(ex,expTab$Experiment[expTab$experiment_id==ex],tr)
for(mod in names(modelNames)) tr <- gsub(mod,modelNames[mod],tr)
tr <- gsub("prior\\\\_","",tr)
for(ab in names(abbv)) tr <- gsub(ab,abbv[ab],tr)

cat(tr,file="results/psModsAppendix.tex")

trOut <- ests1%>%map(~.$outMod)%>%texreg(single.row=TRUE,stars=numeric(0),caption="Coefficient estimates for the outcome regressions described in Section \\ref{sec:application}",label="tab:outModAppendix") ## put in texreg options
for(ex in 1:3) trOut <- gsub(expTab$experiment_id[ex],expTab$Experiment[ex],trOut)
trOut <- gsub("Sp","R",trOut)
for(mod in names(modelNames)) trOut <- gsub(mod,modelNames[mod],trOut)
trOut <- gsub("prior\\\\_","",trOut)
for(ab in names(abbv)) trOut <- gsub(ab,abbv[ab],trOut)

cat(trOut,file="results/outModsAppendix.tex")

