library(readr)
library(magrittr)
library(dplyr)

# Load data
# Load data
data<-readRDS("KellerData/DEAP-data-download_ThompsonScores_ResFamIncome.rds")

F1_res <- read.csv("Residual_ABCD_EUR_Mallard_Factor_1_PRS.csv") #5815
F2_res <- read.csv("Residual_ABCD_EUR_Mallard_Factor_2_PRS.csv") #5815

F1_F2<-merge(F1_res,F2_res,by="src_subject_id") #merge F1 and F2 together
F1_F2$src_subject_id <- gsub(pattern="INV",replacement="NDAR_INV", F1_F2$src_subject_id) #change ID format

#for simplicity, rename PRS column titles
colnames(F1_F2)<-gsub("Residual_PRS_1","PRS_1",colnames(F1_F2))
colnames(F1_F2)<-gsub("Residual_PRS_2","PRS_2",colnames(F1_F2))

#merge in F1 and F2
data.PRS <- merge(data,F1_F2[, c("src_subject_id",'PRS_1','PRS_2')],by="src_subject_id")

#Only keep baseline data
newdata.baseline<-data.PRS[data.PRS$event_name=="baseline_year_1_arm_1",] #5815
colnames(newdata.baseline)[1]<-"subjectkey"

# Load in the spatial extent of each PFN (they're in order, 1-17)
pfn_sizes <- read.csv("KellerData/All_PFN_sizes.csv") #N=9132
# Change PFN column names so they're sensible
colnames(pfn_sizes)<-c("subjectkey","PFN1","PFN2","PFN3","PFN4","PFN5","PFN6","PFN7","PFN8","PFN9","PFN10","PFN11","PFN12","PFN13","PFN14","PFN15","PFN16","PFN17")

#merge PFN sizes with data
all.data <- merge(pfn_sizes,newdata.baseline[,c("subjectkey", setdiff(colnames(newdata.baseline),colnames(pfn_sizes)))],by="subjectkey") #N=4727

# Remove subjects based on the following criteria:
## 1) 8min of retained TRs (600 TRs)
## 2) ABCD Booleans for rest and T1
num_TRs <- read.csv("KellerData/num_TRs_by_subj_redo.csv")
data<-merge(all.data,num_TRs,by="subjectkey")
data.clean<-data[data$numTRs>=600,] #N=4379

# Remove participants based on booleans from ABCD
abcd_imgincl01 <- read.csv("KellerData/abcd_imgincl01.csv")
abcd_imgincl01 <- abcd_imgincl01[abcd_imgincl01$eventnam=="baseline_year_1_arm_1",]
abcd_imgincl01 <- abcd_imgincl01[!abcd_imgincl01$visit=="",] 
abcd_imgincl01 <- abcd_imgincl01[abcd_imgincl01$imgincl_t1w_include==1,] 
abcd_imgincl01 <- abcd_imgincl01[abcd_imgincl01$imgincl_rsfmri_include==1,]
combined.data <- merge(data.clean,abcd_imgincl01[, c("subjectkey",'imgincl_t1w_include')],by="subjectkey") #N=3990

# Add in Family and Site covariates
# Remember to add in this variable about family structure so we can use it as a covariate
family <-read.table("KellerData/acspsw03.txt",header=TRUE)
family.baseline<-family[family$eventname=="baseline_year_1_arm_1",]
abcd.data.almost <- merge(combined.data,family.baseline[, c("subjectkey", setdiff(colnames(family.baseline),colnames(combined.data)))],by="subjectkey")
#N=3990

# also add in the variable for site ID to use as a covariate
site_info <- readRDS("KellerData/DEAP-siteID.rds")
site.baseline<-site_info[site_info$event_name=="baseline_year_1_arm_1",]
colnames(site.baseline)[1]<-"subjectkey"
abcd.data.almost2 <- merge(abcd.data.almost,site.baseline[,c("subjectkey", setdiff(colnames(site.baseline),colnames(abcd.data.almost)))],by="subjectkey")
#N=3990

# Add mean FD
meanFD <- read.csv("KellerData/meanFD_031822.csv")
colnames(meanFD)<-c("subjectkey","meanFD")
meanFD$subjectkey <- gsub(pattern="sub-NDAR",replacement="NDAR_", meanFD$subjectkey)
abcd.data <- merge(abcd.data.almost2,meanFD,by="subjectkey")
#N=3984

# Check for duplicates
abcd.data<-abcd.data[!duplicated(abcd.data$subjectkey),] #N=3982

# Load train/test split from UMinn
traintest<-read_tsv("KellerData/participants.tsv")
traintest.baseline<-traintest[traintest$session_id=="ses-baselineYear1Arm1",c("participant_id","matched_group")]
colnames(traintest.baseline)[1]<-c("subjectkey")
traintest.baseline$subjectkey <- gsub(pattern="sub-NDAR",replacement="NDAR_", traintest.baseline$subjectkey)
traintest.baseline<-traintest.baseline %>% distinct()
abcd.data.traintest <- merge(abcd.data,traintest.baseline,by="subjectkey")
abcd.data.train <- abcd.data.traintest[abcd.data.traintest$matched_group==1,]
abcd.data.test <- abcd.data.traintest[abcd.data.traintest$matched_group==2,]
abcd.data.for.ridge<-abcd.data.traintest[,c("subjectkey","PRS_1","PRS_2","matched_group","interview_age","sex","meanFD","abcd_site","rel_family_id")]
abcd.data.for.ridge.complete<-abcd.data.for.ridge[complete.cases(abcd.data.for.ridge),] #N=3982

abcd.data.for.ridge.complete$subjectkey <- gsub(pattern="NDAR_",replacement="sub-NDAR", abcd.data.for.ridge.complete$subjectkey)
write.csv(abcd.data.for.ridge.complete,"data_for_ridge_PRS.csv")


