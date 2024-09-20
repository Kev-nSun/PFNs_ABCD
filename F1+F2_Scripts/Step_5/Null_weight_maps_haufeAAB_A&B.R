
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "Weight_Haufe_Nulls_041024/";
dirout <- "Null_weight_maps/"

PFNs_hardparcel<-read_cifti('hardparcel_group.dscalar.nii')

for (catg in c('PRS_1','PRS_2')) { 
  count <- 0
  for (dir in c('A','B','C','D','E','F','G','H','I','J')) { #cycle through 10 directories
    for (fold in 0:99) { #cycle through 100 files
      count <- count + 1 #count to determine name of file (1-1000)
      for (grp in c('A','B')) {
      
        infile <- paste0(dirin, 'Weight_Haufe_Nulls_041024_All_cov_', dir, '/', catg, '_', grp, '_Weight_haufetrans_all_', fold,'.mat')
        outfile <- paste0("./",dirout,catg, '_', grp, '_Weight_haufetrans_all_med_wall_', count)
        
        # read in weights_by_vertex.mat
        readin <- readMat(infile);
        read_weights <- readin$Weights
        
        #absolute value of weights
        abs_weights <- abs(read_weights)  
        
        #sum abs weights across 17 networks
        summed_weights <- matrix(0,59412)
        for (i in c(1:59412)) {
          for (f in c(1:17)) {
            summed_weights[i] <- summed_weights[i]+abs_weights[i+(f-1)*59412]
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
      
        if (grp == 'A') {
          vertex_weights_lh_A <- matrix(0,29696,1)
          vertex_weights_lh_A <- vertex_weights_lh
          vertex_weights_rh_A <- matrix(0,29716,1)
          vertex_weights_rh_A <- vertex_weights_rh 
        }
        if (grp == 'B') {
          vertex_weights_lh_B <- matrix(0,29696,1)
          vertex_weights_lh_B <- vertex_weights_lh
          vertex_weights_rh_B <- matrix(0,29716,1)
          vertex_weights_rh_B <- vertex_weights_rh 
        }
      
        #define map variable using hard parcel cifti file
        Weight_map <- PFNs_hardparcel
      
        #assign weight values into hard parcel
        Weight_map$data$cortex_left <- vertex_weights_lh
        Weight_map$data$cortex_right <- vertex_weights_rh
        
        #add back in medial wall
        Weight_map_med <- move_from_mwall(Weight_map)
      
        #save out cifti files in A or B dir
        write_cifti(Weight_map_med,outfile)
      }
      print(count)
    }
  }
}