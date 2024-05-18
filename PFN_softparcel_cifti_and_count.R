
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

PFNs_hardparcel<-read_cifti('hardparcel_group.dscalar.nii') #hard parcel cifti used as template
readin<-readMat('initResamp.mat') #soft parcel mat file
softparcel<-readin$initV

PFN_sizes <- matrix(0,17)

for (net in c(1:17)) {
  
  #assign network values to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_weights_lh[i,1] <- softparcel[i,net]
  }
  
  #assign network values to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_weights_rh[i,1] <- softparcel[i+29696,net]
  }

  #define map variable using hard parcel cifti file
  Soft_map <- PFNs_hardparcel

  #assign network values for each hemi into cifti 
  Soft_map$data$cortex_left <- vertex_weights_lh
  Soft_map$data$cortex_right <- vertex_weights_rh
  
  outfile <- paste0('PFN', net, '_soft_parcel')
  write_cifti(Soft_map,outfile) # save out cifti for each PFN
  
  PFN_count <- 0
  for (vert in c(1:59412)) { #count # of non-0 values in each network soft parcel map
    if (softparcel[vert,net] != 0) {
      PFN_count <- PFN_count + 1
    }
  }
  PFN_sizes[net] <- PFN_count #add PFN_count into matrix of 17
}

write.csv(PFN_sizes, file = "PFN_vert_count.csv") #save out values