library(tidyverse)

dat <- read_csv("data/processed_experiment_data.csv")

as <- dat%>%filter(independent_measure_pair=="Answer Only vs. Scaffolding")

raw <- map_dfr(
  unique(as$experiment_id),
  ~read_csv(
    paste0(
      "~/OneDrive - Worcester Polytechnic Institute (wpi.edu)/iesJohannReloop/ethenExperiments/covariates/",
      .x,"/exp_plogs.csv"))%>%#,col_types = "???????c??????")%>%
    mutate(exp_id=.x))

datraw <- inner_join(as,raw)
with(datraw,table(condition,problem_condition))
with(datraw,table(condition,problem_condition,experiment_id))
with(datraw,table(condition,scaffold_problems_available))
with(datraw,table(condition,scaffold_problems_given))

rawStud <- datraw%>%
  #filter(problem_condition!="Unknown")%>%
  group_by(experiment_id,student_id)%>%
  summarize(
      across(
          c(
              starts_with("scaff"),
              starts_with("hint"),
              starts_with("explanation"),
              answer_given),
          ~sum(.)),
    nrow=n(),
    nprobPart=n_distinct(problem_id,problem_part),
    nprob=n_distinct(problem_id))

dat <- left_join(as,rawStud)

dat$S <- dat$scaffold_problems_given>0
dat$S2 <- dat$scaffold_problems_given>2

dat <- dat%>%group_by(experiment_id)%>%mutate(Srate=mean(S))%>%
  filter(Srate>0)%>%ungroup()%>%select(-Srate)

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


write_csv(dat,"data/as.csv")
