library(readr)
library(ggtext)
library(dplyr)

# First load in the spatial extent of each PFN (they're in order, 1-17)
# Then merge with the same datasets as above (NIH toolbox, trauma, etc.)
pfn_sizes <- read.csv("../KellerData/All_PFN_sizes.csv") #N=9132
# Change PFN column names so they're sensible
colnames(pfn_sizes)<-c("subjectkey","PFN1","PFN2","PFN3","PFN4","PFN5","PFN6","PFN7","PFN8","PFN9","PFN10","PFN11","PFN12","PFN13","PFN14","PFN15","PFN16","PFN17")

# Load data (BASELINE BIFACTOR)
newdata<-readRDS("DEAP-data-download_ABCD_psychopathology_bifactor_scores_BASELINE.rds")
newdata.baseline<-newdata[newdata$event_name=="baseline_year_1_arm_1",]
colnames(newdata.baseline)[1]<-"subjectkey"
all.data <- merge(pfn_sizes,newdata.baseline[,c("subjectkey", setdiff(colnames(newdata.baseline),colnames(pfn_sizes)))],by="subjectkey") #N=9131


# Remove subjects based on the following criteria:
## 1) 8min of retained TRs (600 TRs)
## 2) ABCD Booleans for rest and T1
num_TRs <- read.csv("../KellerData/num_TRs_by_subj_redo.csv")
data<-merge(all.data,num_TRs,by="subjectkey")
data.clean<-data[data$numTRs>=600,] #N=8292

# Remove participants based on booleans from ABCD
abcd_imgincl01 <- read.csv("../KellerData/abcd_imgincl01.csv")
abcd_imgincl01 <- abcd_imgincl01[abcd_imgincl01$eventnam=="baseline_year_1_arm_1",]
abcd_imgincl01 <- abcd_imgincl01[!abcd_imgincl01$visit=="",] 
abcd_imgincl01 <- abcd_imgincl01[abcd_imgincl01$imgincl_t1w_include==1,] 
abcd_imgincl01 <- abcd_imgincl01[abcd_imgincl01$imgincl_rsfmri_include==1,]
combined.data <- merge(data.clean,abcd_imgincl01[, c("subjectkey",'imgincl_t1w_include')],by="subjectkey") #N=7474

# for simplicity, rename general P
colnames(combined.data)<-gsub("General_p","General_PB1",colnames(combined.data))


# Add in Family and Site covariates
# Remember to add in this variable about family structure so we can use it as a covariate
family <-read.table("../KellerData/acspsw03.txt",header=TRUE)
family.baseline<-family[family$eventname=="baseline_year_1_arm_1",]
abcd.data.almost <- merge(combined.data,family.baseline[, c("subjectkey", setdiff(colnames(family.baseline),colnames(combined.data)))],by="subjectkey")
#N=7474

# also add in the variable for site ID to use as a covariate
site_info <- readRDS("../KellerData/DEAP-siteID.rds")
site.baseline<-site_info[site_info$event_name=="baseline_year_1_arm_1",]
colnames(site.baseline)[1]<-"subjectkey"
abcd.data.almost2 <- merge(abcd.data.almost,site.baseline[,c("subjectkey", setdiff(colnames(site.baseline),colnames(abcd.data.almost)))],by="subjectkey")
#N=7474

# Add mean FD
meanFD <- read.csv("../KellerData/meanFD_031822.csv")
colnames(meanFD)<-c("subjectkey","meanFD")
meanFD$subjectkey <- gsub(pattern="sub-NDAR",replacement="NDAR_", meanFD$subjectkey)
abcd.data <- merge(abcd.data.almost2,meanFD,by="subjectkey")
#N=7461

# Check for duplicates
abcd.data<-abcd.data[!duplicated(abcd.data$subjectkey),] # 7459

# Load train/test split from UMinn
traintest<-read_tsv("../KellerData/participants.tsv")
traintest.baseline<-traintest[traintest$session_id=="ses-baselineYear1Arm1",c("participant_id","matched_group")]
colnames(traintest.baseline)[1]<-c("subjectkey")
traintest.baseline$subjectkey <- gsub(pattern="sub-NDAR",replacement="NDAR_", traintest.baseline$subjectkey)
traintest.baseline<-traintest.baseline %>% distinct()
abcd.data.traintest <- merge(abcd.data,traintest.baseline,by="subjectkey")
abcd.data.train <- abcd.data.traintest[abcd.data.traintest$matched_group==1,]
abcd.data.test <- abcd.data.traintest[abcd.data.traintest$matched_group==2,]
abcd.data.for.ridge<-abcd.data.traintest[,c("subjectkey","General_PB1","matched_group","interview_age","sex","meanFD","abcd_site","rel_family_id")]
abcd.data.for.ridge.complete<-abcd.data.for.ridge[complete.cases(abcd.data.for.ridge),]

abcd.data.for.ridge.complete$subjectkey <- gsub(pattern="NDAR_",replacement="sub-NDAR", abcd.data.for.ridge.complete$subjectkey)
write.csv(abcd.data.for.ridge.complete,"data_for_ridge_pfac_BASELINE.csv")
