
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "../../Step_5/BASELINE/Weight_Haufe_031324/";
dirout <- "Rel_max_maps/"
dir.create(dirout)

PFNs_hardparcel<-read_cifti('../hardparcel_group.dscalar.nii')

#read in A and B weights separately
catg <- 'General_PB1'
for (grp in c('A','B')) {
  
  infile <- paste0(dirin, catg, '_', grp, '_Weight_haufetrans_all','.mat')
  
  # read in weights_by_vertex.mat
  readin <- readMat(infile);
  read_weights <- readin$PB.test
  
  #format Haufe weights into a matrix of 17 PFNs
  all_weights <- matrix(0,59412,17)
  for (vert in c(1:59412)) {
    for (PFN in c(1:17)) {
      all_weights[vert,PFN] <- read_weights[vert+(PFN-1)*59412]
    }
  }
  
  if (grp == 'A') {
    all_weights_A <- matrix(0,59412,17)
    all_weights_A <- all_weights
    concat_weights_A <- read_weights
  }
  if (grp == 'B') {
    all_weights_B <- matrix(0,59412,17)
    all_weights_B <- all_weights
    concat_weights_B <- read_weights
  }
}

#average Haufe weights between groups A and B
all_weights_avg <- matrix(0,59412,17) 
for (vert in c(1:59412)) {
  for (PFN in c(1:17)) {
    all_weights_avg[vert,PFN] <- sum(all_weights_A[vert,PFN]+all_weights_B[vert,PFN])/2
  }
}

#Max Pos/Neg Weight Maps
max_pos_weight <- matrix(0,59412,2)
max_neg_weight <- matrix(0,59412,2)

#determine max pos and neg weight and PFN ID for each vertex
for (vert in c(1:59412)) {
  max_tracker_pos <- 0
  max_tracker_neg <- 0
  for (PFN in c(1:17)) {
    if (all_weights_avg[vert,PFN] > max_tracker_pos) { #check if PFN weight is larger than max pos weight tracker
      max_tracker_pos <- all_weights_avg[vert,PFN]
      pos_PFNID <- PFN
    } else if (all_weights_avg[vert,PFN] < max_tracker_neg) { #check if PFN weight is smaller than max neg weight tracker
      max_tracker_neg <- all_weights_avg[vert,PFN]
      neg_PFNID <- PFN
    }
  }
  max_pos_weight[vert,1] <- max_tracker_pos
  max_pos_weight[vert,2] <- pos_PFNID
  max_neg_weight[vert,1] <- max_tracker_neg
  max_neg_weight[vert,2] <- neg_PFNID
}
write.csv(max_pos_weight,file = paste0(dirout,catg,"_max_pos_weights.csv"))
write.csv(max_neg_weight,file = paste0(dirout,catg,"_max_neg_weights.csv"))

#determine relative threshold for max pos and max neg weights (top 1% by magnitude)
pos_thres <- quantile(max_pos_weight[,1],prob = 0.95) #10% = 0.00374, 5% = 0.00455, 1% = 0.00675
neg_thres <- quantile(max_neg_weight[,1],prob = 0.05) #10% = -0.00303, 5% = -0.00376, 1% = -0.00541


#MAX POS MAP
#assign IDs to left hemi
vertex_weight_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weight_lh[i,1] <- max_pos_weight[i,1]
}
#assign IDs to right hemi
vertex_weight_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weight_rh[i,1] <- max_pos_weight[i+29696,1]
}
#Save out pre-clustered cifti
Pos_weight_map <- PFNs_hardparcel
#assign weight values into hard parcel
Pos_weight_map$data$cortex_left <- vertex_weight_lh
Pos_weight_map$data$cortex_right <- vertex_weight_rh
#define output file for averaged map
outfile <- paste0("./",dirout,catg,'_preclus_max_pos')
#save out cifti files
write_cifti(Pos_weight_map,outfile)

#find clusters using pos_thres and 25 mm2 surface area in cifti
infile <- outfile
clust_file <- paste0(dirout,catg,'_cluster_pos_5.dscalar.nii')
cmd <- paste0('wb_command -cifti-find-clusters ',infile,'.dscalar.nii ',pos_thres,' ',25,' ',0,' ',0,' COLUMN ',clust_file,' -left-surface ../../Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii -right-surface ../../Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii')
system(cmd)

#read in clustered vertices and relabel them with PFNID
clus_max_pos<-read_cifti(clust_file)
#L hemi
vertex_ID_lh <- matrix(0,29696)
clusts_lh <- clus_max_pos$data$cortex_left
for (i in c(1:29696)) {
  if (clusts_lh[i] > 0){ #if vertex passed value + cluster threshold
    vertex_ID_lh[i] <- max_pos_weight[i,2] #assign vertex it's respective PFN ID
  }
}
#R hemi
vertex_ID_rh <- matrix(0,29716)
clusts_rh <- clus_max_pos$data$cortex_right
for (i in c(1:29716)) {
  if (clusts_rh[i] > 0){ #if vertex passed value + cluster threshold
    vertex_ID_rh[i] <- max_pos_weight[i+29696,2]
  }
}

#Save clustered left + right hem to mat files, converted to gifit--> cifti in matlab "Labeling_PFNs_into_cifti.m"
# For L hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_lh <- vertex_ID_lh
writeMat(paste0(dirout,"Max_Pos_L_clust_5.mat"), x = names_list)
# For R hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_rh <- vertex_ID_rh
writeMat(paste0(dirout,"Max_Pos_R_clust_5.mat"), x = names_list)



#MAX NEG MAP
#assign IDs to left hemi
vertex_weight_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weight_lh[i,1] <- max_neg_weight[i,1]
}
#assign IDs to right hemi
vertex_weight_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weight_rh[i,1] <- max_neg_weight[i+29696,1]
}
#Save out pre-clustered cifti
Neg_weight_map <- PFNs_hardparcel
#assign weight values into hard parcel
Neg_weight_map$data$cortex_left <- vertex_weight_lh
Neg_weight_map$data$cortex_right <- vertex_weight_rh
#define output file for averaged map
outfile <- paste0("./",dirout,catg,'_preclus_max_neg')
#save out cifti files
write_cifti(Neg_weight_map,outfile)

#find clusters using neg_thres and 25 mm2 surface area in cifti
infile <- outfile
clust_file <- paste0(dirout,catg,'_cluster_neg_5.dscalar.nii')
cmd <- paste0('wb_command -cifti-find-clusters ',infile,'.dscalar.nii ',neg_thres,' ',25,' ',0,' ',0,' COLUMN ',clust_file,' -less-than -left-surface ../../Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii -right-surface ../../Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii')
system(cmd)

#read in clustered vertices and relabel them with PFNID
clus_max_neg<-read_cifti(clust_file)
#L hemi
vertex_ID_lh <- matrix(0,29696)
clusts_lh <- clus_max_neg$data$cortex_left
for (i in c(1:29696)) {
  if (clusts_lh[i] > 0){ #if vertex passed value + cluster threshold
    vertex_ID_lh[i] <- max_neg_weight[i,2] #assign vertex it's respective PFN ID
  }
}
#R hemi
vertex_ID_rh <- matrix(0,29716)
clusts_rh <- clus_max_neg$data$cortex_right
for (i in c(1:29716)) {
  if (clusts_rh[i] > 0){ #if vertex passed value + cluster threshold
    vertex_ID_rh[i] <- max_neg_weight[i+29696,2]
  }
}

#Save clustered left + right hem to mat files, converted to gifit--> cifti in matlab "Labeling_PFNs_into_cifti.m"
# For L hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_lh <- vertex_ID_lh
writeMat(paste0(dirout,"Max_Neg_L_clust_5.mat"), x = names_list)
# For R hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_rh <- vertex_ID_rh
writeMat(paste0(dirout,"Max_Neg_R_clust_5.mat"), x = names_list)