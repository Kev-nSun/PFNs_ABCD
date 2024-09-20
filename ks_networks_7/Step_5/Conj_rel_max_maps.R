# Conjunction maps between max pos and neg p-fac, F1, and F2

library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

#read in rel max maps
Pfac_pos <- read_cifti('Rel_max_maps/Final/Pfac_Max_Pos_clust_10_mwall.dlabel.nii')
Pfac_neg <- read_cifti('Rel_max_maps/Final/Pfac_Max_Neg_clust_10_mwall.dlabel.nii')
F1_pos <- read_cifti('Rel_max_maps/Final/Max_Pos_PRS_1_clust_10_mwall.dlabel.nii')
F1_neg <- read_cifti('Rel_max_maps/Final/Max_Neg_PRS_1_clust_10_mwall.dlabel.nii')
F2_pos <- read_cifti('Rel_max_maps/Final/Max_Pos_PRS_2_clust_10_mwall.dlabel.nii')
F2_neg <- read_cifti('Rel_max_maps/Final/Max_Neg_PRS_2_clust_10_mwall.dlabel.nii')

#PFAC AND F1
#loop through L and R hemis, comparing pfac_pos to F1_pos PFN IDs
Pfac_F1_pos_map <- Pfac_pos

Pfac_F1_pos_counter <- 0
Pfac_F1_pos_L <- matrix(0,29696)
for (vert in c(1:29696)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_pos$data$cortex_left[vert] != 0 & F1_pos$data$cortex_left[vert] != 0 & Pfac_pos$data$cortex_left[vert] == F1_pos$data$cortex_left[vert]) {
    Pfac_F1_pos_L[vert] <- Pfac_pos$data$cortex_left[vert]
    Pfac_F1_pos_counter <- Pfac_F1_pos_counter + 1
  }
}
Pfac_F1_pos_map$data$cortex_left <- Pfac_F1_pos_L

Pfac_F1_pos_R <- matrix(0,29716)
for (vert in c(1:29716)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_pos$data$cortex_right[vert] != 0 & F1_pos$data$cortex_right[vert] != 0 & Pfac_pos$data$cortex_right[vert] == F1_pos$data$cortex_right[vert]) {
    Pfac_F1_pos_R[vert] <- Pfac_pos$data$cortex_right[vert]
    Pfac_F1_pos_counter <- Pfac_F1_pos_counter + 1
  }
}
Pfac_F1_pos_map$data$cortex_right <- Pfac_F1_pos_R

#save out new cifti file
dir <- 'Rel_max_maps/Final/'
catg<- 'pfac_F1'
sign<- "pos"
outfile <- paste0(dir,'Conj_',catg,'_',sign,'_clust_10')
#write_cifti(Pfac_F1_pos_map,outfile)


#loop through L and R hemis, comparing pfac_neg to F1_neg PFN IDs
Pfac_F1_neg_map <- Pfac_neg

Pfac_F1_neg_counter <- 0
Pfac_F1_neg_L <- matrix(0,29696)
for (vert in c(1:29696)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_neg$data$cortex_left[vert] != 0 & F1_neg$data$cortex_left[vert] != 0 & Pfac_neg$data$cortex_left[vert] == F1_neg$data$cortex_left[vert]) {
    Pfac_F1_neg_L[vert] <- Pfac_neg$data$cortex_left[vert]
    Pfac_F1_neg_counter <- Pfac_F1_neg_counter + 1
  }
}
Pfac_F1_neg_map$data$cortex_left <- Pfac_F1_neg_L

Pfac_F1_neg_R <- matrix(0,29716)
for (vert in c(1:29716)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_neg$data$cortex_right[vert] != 0 & F1_neg$data$cortex_right[vert] != 0 & Pfac_neg$data$cortex_right[vert] == F1_neg$data$cortex_right[vert]) {
    Pfac_F1_neg_R[vert] <- Pfac_neg$data$cortex_right[vert]
    Pfac_F1_neg_counter <- Pfac_F1_neg_counter + 1
  }
}
Pfac_F1_neg_map$data$cortex_right <- Pfac_F1_neg_R

#save out new cifti file
dir <- 'Rel_max_maps/Final/'
catg<- 'pfac_F1'
sign<- "neg"
outfile <- paste0(dir,'Conj_',catg,'_',sign,'_clust_10')
#write_cifti(Pfac_F1_neg_map,outfile)



#PFAC AND F2
#loop through L and R hemis, comparing pfac_pos to F2_pos PFN IDs
Pfac_F2_pos_map <- Pfac_pos

Pfac_F2_pos_counter <- 0
Pfac_F2_pos_L <- matrix(0,29696)
for (vert in c(1:29696)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_pos$data$cortex_left[vert] != 0 & F2_pos$data$cortex_left[vert] != 0 & Pfac_pos$data$cortex_left[vert] == F2_pos$data$cortex_left[vert]) {
    Pfac_F2_pos_L[vert] <- Pfac_pos$data$cortex_left[vert]
    Pfac_F2_pos_counter <- Pfac_F2_pos_counter + 1
  }
}
Pfac_F2_pos_map$data$cortex_left <- Pfac_F2_pos_L

Pfac_F2_pos_R <- matrix(0,29716)
for (vert in c(1:29716)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_pos$data$cortex_right[vert] != 0 & F2_pos$data$cortex_right[vert] != 0 & Pfac_pos$data$cortex_right[vert] == F2_pos$data$cortex_right[vert]) {
    Pfac_F2_pos_R[vert] <- Pfac_pos$data$cortex_right[vert]
    Pfac_F2_pos_counter <- Pfac_F2_pos_counter + 1
  }
}
Pfac_F2_pos_map$data$cortex_right <- Pfac_F2_pos_R

#save out new cifti file
dir <- 'Rel_max_maps/Final/'
catg<- 'pfac_F2'
sign<- "pos"
outfile <- paste0(dir,'Conj_',catg,'_',sign,'_clust_10')
#write_cifti(Pfac_F2_pos_map,outfile)


#loop through L and R hemis, comparing pfac_neg to F2_neg PFN IDs
Pfac_F2_neg_map <- Pfac_neg

Pfac_F2_neg_counter <- 0
Pfac_F2_neg_L <- matrix(0,29696)
for (vert in c(1:29696)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_neg$data$cortex_left[vert] != 0 & F2_neg$data$cortex_left[vert] != 0 & Pfac_neg$data$cortex_left[vert] == F2_neg$data$cortex_left[vert]) {
    Pfac_F2_neg_L[vert] <- Pfac_neg$data$cortex_left[vert]
    Pfac_F2_neg_counter <- Pfac_F2_neg_counter + 1
  }
}
Pfac_F2_neg_map$data$cortex_left <- Pfac_F2_neg_L

Pfac_F2_neg_R <- matrix(0,29716)
for (vert in c(1:29716)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (Pfac_neg$data$cortex_right[vert] != 0 & F2_neg$data$cortex_right[vert] != 0 & Pfac_neg$data$cortex_right[vert] == F2_neg$data$cortex_right[vert]) {
    Pfac_F2_neg_R[vert] <- Pfac_neg$data$cortex_right[vert]
    Pfac_F2_neg_counter <- Pfac_F2_neg_counter + 1
  }
}
Pfac_F2_neg_map$data$cortex_right <- Pfac_F2_neg_R

#save out new cifti file
dir <- 'Rel_max_maps/Final/'
catg<- 'pfac_F2'
sign<- "neg"
outfile <- paste0(dir,'Conj_',catg,'_',sign,'_clust_10')
#write_cifti(Pfac_F2_neg_map,outfile)



#F1 AND F2
#loop through L and R hemis, comparing F1_pos to F2_pos PFN IDs
F1_F2_pos_map <- F1_pos

F1_F2_pos_counter <- 0
F1_F2_pos_L <- matrix(0,29696)
for (vert in c(1:29696)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (F1_pos$data$cortex_left[vert] != 0 & F2_pos$data$cortex_left[vert] != 0 & F1_pos$data$cortex_left[vert] == F2_pos$data$cortex_left[vert]) {
    F1_F2_pos_L[vert] <- F1_pos$data$cortex_left[vert]
    F1_F2_pos_counter <- F1_F2_pos_counter + 1
  }
}
F1_F2_pos_map$data$cortex_left <- F1_F2_pos_L

F1_F2_pos_R <- matrix(0,29716)
for (vert in c(1:29716)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (F1_pos$data$cortex_right[vert] != 0 & F2_pos$data$cortex_right[vert] != 0 & F1_pos$data$cortex_right[vert] == F2_pos$data$cortex_right[vert]) {
    F1_F2_pos_R[vert] <- F1_pos$data$cortex_right[vert]
    F1_F2_pos_counter <- F1_F2_pos_counter + 1
  }
}
F1_F2_pos_map$data$cortex_right <- F1_F2_pos_R

#save out new cifti file
dir <- 'Rel_max_maps/Final/'
catg<- 'F1_F2'
sign<- "pos"
outfile <- paste0(dir,'Conj_',catg,'_',sign,'_clust_10')
#write_cifti(F1_F2_pos_map,outfile)


#loop through L and R hemis, comparing F1_neg to F2_neg PFN IDs
F1_F2_neg_map <- F1_neg

F1_F2_neg_counter <- 0
F1_F2_neg_L <- matrix(0,29696)
for (vert in c(1:29696)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (F1_neg$data$cortex_left[vert] != 0 & F2_neg$data$cortex_left[vert] != 0 & F1_neg$data$cortex_left[vert] == F2_neg$data$cortex_left[vert]) {
    F1_F2_neg_L[vert] <- F1_neg$data$cortex_left[vert]
    F1_F2_neg_counter <- F1_F2_neg_counter + 1
  }
}
F1_F2_neg_map$data$cortex_left <- F1_F2_neg_L

F1_F2_neg_R <- matrix(0,29716)
for (vert in c(1:29716)) { #if vertex value !=0 and two map PFNIDs are equal, assign PFNID to new map
  if (F1_neg$data$cortex_right[vert] != 0 & F2_neg$data$cortex_right[vert] != 0 & F1_neg$data$cortex_right[vert] == F2_neg$data$cortex_right[vert]) {
    F1_F2_neg_R[vert] <- F1_neg$data$cortex_right[vert]
    F1_F2_neg_counter <- F1_F2_neg_counter + 1
  }
}
F1_F2_neg_map$data$cortex_right <- F1_F2_neg_R

#save out new cifti file
dir <- 'Rel_max_maps/Final/'
catg<- 'F1_F2'
sign<- "neg"
outfile <- paste0(dir,'Conj_',catg,'_',sign,'_clust_10')
#write_cifti(F1_F2_neg_map,outfile)

pos_overlap_count <- c(Pfac_F1_pos_counter,Pfac_F2_pos_counter,F1_F2_pos_counter)
neg_overlap_count <- c(Pfac_F1_neg_counter,Pfac_F2_neg_counter,F1_F2_neg_counter)
overlap_count <- matrix(0,2,3)
colnames(overlap_count) <- c("Pfac_F1","Pfac_F2","F1_F2")
rownames(overlap_count) <- c("Pos_overlap","Neg_overlap")
overlap_count[1,] <- pos_overlap_count
overlap_count[2,] <- neg_overlap_count
write.csv(overlap_count, file = "Conj_maps_overlap_count.csv")


