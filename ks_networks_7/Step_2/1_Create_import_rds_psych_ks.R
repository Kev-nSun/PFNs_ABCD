library(readr)
library(magrittr)
library(dplyr)

# Load data
data<-readRDS("../KellerData/DEAP-data-download_ThompsonScores_ResFamIncome.rds")

psyc <- read.csv("ABCD_psychopathology_bifactor_scores_BASELINE_NDAR.csv")

colnames(psyc)[1] <- "src_subject_id"
data.psyc <- merge(data,psyc[, c("src_subject_id",'Factor1','Factor2','Factor3','Factor4','Factor5','Factor6','Factor7','Factor8','General_p')],by="src_subject_id")

write.csv(data.psyc,"DEAP-data-download_ABCD_psychopathology_bifactor_scores_BASELINE.csv")
saveRDS(data.psyc, file = "DEAP-data-download_ABCD_psychopathology_bifactor_scores_BASELINE.rds")
