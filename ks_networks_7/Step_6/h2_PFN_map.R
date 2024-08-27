
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)
library(dplyr)

dirout <- "h2_map/"
dir.create(dirout)

PFNs_hardparcel<-read_cifti('../Step_6/hardparcel_group.dscalar.nii')
PFN_h2_non_zero<-read.csv('PFN_7459_ETD_UNI_age_sex_3.19.labeled.csv')

h2_summary<-select(PFN_h2_non_zero,var,a2,pNoA,pNoA_fdr,log10_pNoA_fdr) #summary of important vals
write.csv(h2_summary, file = "h2_summary.csv") #save out summary

PFN_h2_all<- matrix(0,1010004,3)
index <- 1
non_zero_vert <- PFN_h2_non_zero[index,3] #PFN vertex id # in non-zero list

for (vert in c(1:1010004)) {
  if (vert == non_zero_vert) { #if PFN vertex id # is == to vert, replace 0's with values for var, a2, and log10_pNoA_fdr
    PFN_h2_all[vert,1]<- PFN_h2_non_zero[index,3]
    PFN_h2_all[vert,2]<- PFN_h2_non_zero[index,11]
    PFN_h2_all[vert,3]<- PFN_h2_non_zero[index,27]
    index <- index + 1 # move to next vertex value in non-zero list if place has been found
    non_zero_vert <- PFN_h2_non_zero[index,3]
  }
  print(vert)
}
#write.csv(PFN_h2_all, file = "All_PFN_h2_pvals.csv") #save out all h2 values and pvals with zeros

#take highest h2 value across 17 networks for each vertex
PFN_h2_top <- matrix(0,59412,2)
for (vert in c(1:59412)) {
  top_h2_tracker = 0 #set h2 tracker to 0
  for (net in c(1:17)) {
    current_h2 <- PFN_h2_all[vert+(net-1)*59412,2] #record current h2
    if (top_h2_tracker < current_h2) #check if current h2 is higher than tracker
    {
      top_h2_tracker = current_h2 #record new highest h2
      PFN_h2_top[vert,1] <- current_h2 #save current highest h2
      PFN_h2_top[vert,2] <- PFN_h2_all[vert+(net-1)*59412,3] #save corresponding log10_pNoA_fdr val
    }
  }
}
#write.csv(PFN_h2_top, file = paste0("./",dirout,"Max_PFN_h2_pvals.csv")) #save out max h2 values with log10_pNoA_fdr values

#Heritability map with ALL vertices
#assign h2 to left hemi
vertex_h2_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_h2_lh[i,1] <- PFN_h2_top[i,1]
}

#assign h2 to right hemi
vertex_h2_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_h2_rh[i,1] <- PFN_h2_top[i+29696,1]
}

#define map variable using hard parcel cifti file
h2_map <- PFNs_hardparcel

#assign weight values into hard parcel
h2_map$data$cortex_left <- vertex_h2_lh
h2_map$data$cortex_right <- vertex_h2_rh

#add back in medial wall
h2_map_med <- move_from_mwall(h2_map)

#define output file for map
outfile <- paste0("./",dirout,'Max_PFN_h2')

#save out cifti files
write_cifti(h2_map_med,outfile)


#Heritability map with ONLY SIG vertices
#assign h2 to left hemi
vertex_h2_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  if (PFN_h2_top[i,2] > 1.30103) { #only add h2 if sig
    vertex_h2_lh[i,1] <- PFN_h2_top[i,1]
  }
}

#assign h2 to right hemi
vertex_h2_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  if (PFN_h2_top[i+29696,2] > 1.30103) { #only add h2 if sig
    vertex_h2_rh[i,1] <- PFN_h2_top[i+29696,1]
  }
}

#define map variable using hard parcel cifti file
h2_map <- PFNs_hardparcel

#assign weight values into hard parcel
h2_map$data$cortex_left <- vertex_h2_lh
h2_map$data$cortex_right <- vertex_h2_rh

#add back in medial wall
h2_map_med <- move_from_mwall(h2_map)

#define output file for map
outfile <- paste0("./",dirout,'Max_PFN_SIG_ONLY_h2')

#save out cifti files
write_cifti(h2_map_med,outfile)


#CALCULATE AVG h2 OF TOP % PFN WEIGHTS OF PRS MODELS
#Check what proportion of top % PFN weights of PRS models have sig h2
dirout <- "h2_top%/"
dir.create(dirout)

PFN_h2_top_perc_read <- read.csv("h2_map/Max_PFN_h2_pvals.csv") #read in max PFN vertex heritabiltiy
PFN_h2_top <- subset(PFN_h2_top_perc_read,select=-c(X))

Avg_h2_top_perc_df <- data.frame() #list of average h2 of top 10% weights for each PRS
sig_h2_vert_counter <- matrix(0,2)
PRS_top_perc_maps <- matrix(0,59412,2)
for (catg in c(1:2)){
  if (catg == 1) {
    Top_perc <- read.csv("./top%_weights/PRS_1_Sum_of_weights_top_5%.csv") #read in top 5% summed PFN weights for each PRS model
  } else {
    Top_perc <- read.csv("./top%_weights/PRS_2_Sum_of_weights_top_5%.csv")
  }
  h2_Top_perc <- data.frame()

  for (vert in c(1:59412)) {
    if (Top_perc[vert,2] > 0) { #if value > 0, vertex is in top 10% of summed weights
      h2_new_val <- PFN_h2_top[vert,1]
      h2_Top_perc <- rbind(h2_Top_perc,h2_new_val) #append h2 value corresponding to list of top 10% for each PRS (to average)
      PRS_top_perc_maps[vert,catg] <- h2_new_val #append h2 value to map of top 10% for each PRS
      if (PFN_h2_top[vert,2] > 1.30103) { #count sig h2 values in top 10% for each PRS
        sig_h2_vert_counter[catg] <- sig_h2_vert_counter[catg] + 1
      }
    }
  }
  colnames(h2_Top_perc) <- c("h2")
  Avg_h2_Top_perc <- mean(h2_Top_perc$h2) #take average of h2 values
  Avg_h2_top_perc_df <- rbind(Avg_h2_top_perc_df,Avg_h2_Top_perc)
}
colnames(Avg_h2_top_perc_df) <- c("Avg_h2")
rownames(Avg_h2_top_perc_df) <- c("PRS_1","PRS_2")
#write.csv(Avg_h2_top_perc_df, file = paste0("./",dirout,"Avg_h2_of_top5_weights.csv")) #save out avg h2 values for top weights of PRS_1 and PRS_2
#write.csv(PRS_top_perc_maps, file = paste0("./",dirout,"h2_of_top5_PRS_weights.csv")) #save out  h2 values for top weights of PRS_1 and PRS_2

#Heritability map with ONLY TOP 5% F1/F2 Vertices
for (catg in c(1:2)){
  #assign h2 to left hemi
  vertex_h2_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_h2_lh[i,1] <- PRS_top_perc_maps[i,catg]
  }
  
  #assign h2 to right hemi
  vertex_h2_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_h2_rh[i,1] <- PRS_top_perc_maps[i+29696,catg]
  }
  
  #define map variable using hard parcel cifti file
  h2_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel
  h2_map$data$cortex_left <- vertex_h2_lh
  h2_map$data$cortex_right <- vertex_h2_rh
  
  #add back in medial wall
  h2_map_med <- move_from_mwall(h2_map)
  
  #define output files for maps
  if (catg == 1) {
    outfile <- paste0("./",dirout,'Max_PFN_h2_F1-top5')
  } else {
    outfile <- paste0("./",dirout,'Max_PFN_h2_F2-top5')
  }
  
  #save out cifti files
  #write_cifti(h2_map_med,outfile)
}


#CALCULATE AVG h2 OF TOP % PFN WEIGHTS OF Pfac
#Check what proportion of top % PFN weights of pfac models have sig h2
dirout <- "h2_top%/"

sig_h2_vert_counter <- 0
Pfac_top_perc_map <- matrix(0,59412)
top_perc <- read.csv("./top%_weights/General_PB1_Sum_of_weights_top_5%.csv") #read in top 5% summed PFN weights for Pfac
h2_top_perc <- data.frame()
  
for (vert in c(1:59412)) {
  if (top_perc[vert,2] > 0) { #if value > 0, vertex is in top 10% of summed weights
    h2_new_val <- PFN_h2_top[vert,1]
    h2_top_perc <- rbind(h2_top_perc,h2_new_val) #append h2 value corresponding to list of top 10% for pfac
    Pfac_top_perc_map[vert] <- h2_new_val #append h2 value to map of top 10% for pfac
    if (PFN_h2_top[vert,2] > 1.30103) { #count sig h2 values in top 10% for RS
      sig_h2_vert_counter <- sig_h2_vert_counter + 1
    }
  }
}
colnames(h2_top_perc) <- c("h2")
Avg_h2_top_perc <- mean(h2_top_perc$h2) #take average of h2 values
write.csv(Avg_h2_top_perc, file = paste0("./",dirout,"Avg_h2_of_top5_pfac_weights.csv")) #save out avg h2 values for top weights of Pfac
write.csv(Pfac_top_perc_map, file = paste0("./",dirout,"h2_of_top5_pfac_weights.csv")) #save out h2 values for top weights of Pfac

#Heritability map with ONLY TOP 5% pfac Vertices
#assign h2 to left hemi
vertex_h2_lh <- matrix(0,29696)
for (i in c(1:29696)) {
  vertex_h2_lh[i] <- Pfac_top_perc_map[i]
}
  
  #assign h2 to right hemi
vertex_h2_rh <- matrix(0,29716)
for (i in c(1:29716)) {
  vertex_h2_rh[i] <- Pfac_top_perc_map[i+29696]
}
  
#define map variable using hard parcel cifti file
h2_map <- PFNs_hardparcel
  
#assign weight values into hard parcel
h2_map$data$cortex_left <- vertex_h2_lh
h2_map$data$cortex_right <- vertex_h2_rh
  
#add back in medial wall
h2_map_med <- move_from_mwall(h2_map)
  
#define output files
outfile <- paste0("./",dirout,'Max_PFN_h2_pfac-top5') 
  
#save out cifti files
write_cifti(h2_map_med,outfile)

