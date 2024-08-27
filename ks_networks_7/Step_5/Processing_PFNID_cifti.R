# Fixes the Matlab output files (*.dlabel.nii) which are missing the medial wall mask values
# Based on hardparcel_group.dscalar.nii as template

library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

#medial wall mask template is hard parcel file
mwall_template <- read_cifti('../hardparcel_group.dscalar.nii')

#read in map that needs medial wall mask values
map_name <- 'Max_Neg'
dir <- 'PFNID_Visualize'
map <- read_cifti(paste0(dir,'/',map_name,'_AtlasLabel.dlabel.nii'))

#add in mwall mask values based on template
map$meta$cortex$medial_wall_mask <- mwall_template$meta$cortex$medial_wall_mask

#save out new cifti file
outfile <- paste0(dir,'/',map_name,'_AtlasLabel_mwall')
write_cifti(map,outfile)


#add medial wall into each PFN map
for ( PFN in c(1:17)) {
  infile <- paste0(paste0(dir,'/',map_name,'_AtlasLabel_Network_',PFN,'.dlabel.nii'));
  outfile <- paste0(paste0(dir,'/',map_name,'_AtlasLabel_Network_',PFN,'_mwall'))
  
  map<-read_cifti(infile);
  
  map$meta$cortex$medial_wall_mask <- mwall_template$meta$cortex$medial_wall_mask
  
  #save out cifti files
  write_cifti(map,outfile)
}