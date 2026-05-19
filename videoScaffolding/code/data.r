
dat0 <- read_csv("rawdata/processed_experiment_data.csv")

as <- dat0%>%filter(independent_measure_pair=="Answer Only vs. Video Scaffolding")

### extract experiment data
## the data downloaded from OSF is a zip file of zip files--this code extracts the needed
## zipped sub-folders from the larger zip file
if(!all(file.exists(paste0("rawdata/",unique(as$experiment_id),".zip"))))
    unzip("rawdata/experiment_dataset_2021-09-23.zip",
          paste0(unique(as$experiment_id),".zip"),
          exdir="rawdata")


raw <- map_dfr(
  unique(as$experiment_id),
  ~read_csv(unz(paste0("rawdata/",.x,".zip"),"exp_plogs.csv"))%>%
  filter(student_id%in%as$student_id[as$experiment_id==.x])
)

raw <- mutate(raw,tut=hints_given+scaffold_problems_given+explanation_given+answer_given)

raw <- filter(raw,experiment_id%in%as$experiment_id,student_id%in%as$student_id)%>%
filter(problem_condition!="Nonexperimental",problem_condition!="Unknown")


datraw <- inner_join(as,raw)

rawStud <- datraw%>%
                                        #filter(problem_condition!="Unknown")%>%
  group_by(experiment_id,student_id,condition)%>%
  summarize(
      nscaff=n_distinct(scaffold_id,na.rm=TRUE),
      nscaffGiven=n_distinct(scaffold_problems_given,na.rm=TRUE),
    nrow=n(),
    nprobPart=n_distinct(problem_id,problem_part),
    nprob=n_distinct(problem_id),
    across(
          c(
              answer_before_tutoring,
              starts_with("scaff"),
              starts_with("hint"),
              starts_with("explanation"),
              answer_given),
          ~sum(.)))

dat <- left_join(as,rawStud)

dat$S <- dat$scaffold_problems_given>0

dat <- dat%>%group_by(experiment_id)%>%mutate(nprobS=scale(nprob)[,1])%>%ungroup()

dat <- dat%>%group_by(experiment_id)%>%mutate(Srate=mean(S))%>%
  filter(Srate>0)%>%ungroup()%>%dplyr::select(-Srate)



write_csv(dat,"data/as.csv")


