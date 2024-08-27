
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "../../Step_5/BASELINE/Haufe_Nulls_031824/";
dirout <- "weight_map_HaufeAAB_nulls/"
dir.create(dirout)

PFNs_hardparcel<-read_cifti('../hardparcel_group.dscalar.nii')

count <- 0
catg <- 'General_PB1'
for (dir in c('A','B','C','D','E','F','G','H','I','J')) { #cycle through 10 directories
  for (fold in 0:99) { #cycle through 100 files
    count <- count + 1 #count to determine placement of values in sum_of_weights null matrices
    for (grp in c('A','B')) {
    
      infile <- paste0(dirin, 'Weight_Haufe_Nulls_031824_', dir, '/', catg, '_', grp, '_Weight_haufetrans_all_', fold,'.mat')
      
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
      }
      if (grp == 'B') {
        all_weights_B <- matrix(0,59412,17)
        all_weights_B <- all_weights
      }
    }
    
    #average Haufe weights between groups A and B
    all_weights_avg <- matrix(0,59412,17) 
    for (vert in c(1:59412)) {
      for (PFN in c(1:17)) {
        all_weights_avg[vert,PFN] <- sum(all_weights_A[vert,PFN]+all_weights_B[vert,PFN])/2
      }
    }
    
    #absolute value of weights
    abs_avg_weights <- abs(all_weights_avg)
    
    #sum abs weights across 17 networks
    summed_weights <- matrix(0,59412)
    for (i in c(1:59412)) {
      for (f in c(1:17)) {
        summed_weights[i] <- summed_weights[i]+abs_avg_weights[i+(f-1)*59412]
      }
    }
    
    #assign summed weights to left hemi
    vertex_weights_lh <- matrix(0,29696,1)
    for (i in c(1:29696)) {
      vertex_weights_lh[i,1] <- summed_weights[i]
    }
    
    #assign summed weights to right hemi
    vertex_weights_rh <- matrix(0,29716,1)
    for (i in c(1:29716)) {
      vertex_weights_rh[i,1] <- summed_weights[i+29696]
    }
    
    #define map variable using hard parcel cifti file
    Weight_map <- PFNs_hardparcel
    
    #assign weight values into hard parcel
    Weight_map$data$cortex_left <- vertex_weights_lh
    Weight_map$data$cortex_right <- vertex_weights_rh
    
    #add back in medial wall
    Weight_map_med <- move_from_mwall(Weight_map)
    
    outfile <- paste0(dirout,catg, '_Avg_Weight_haufetrans_all_med_wall_', count)
    
    #save out cifti files
    write_cifti(Weight_map_med,outfile)
    print(count)
  }
}