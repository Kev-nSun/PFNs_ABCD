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
map_name <- 'Max_Neg_PRS_2_clust_5'
dir <- 'Rel_max_maps'
map <- read_cifti(paste0(dir,'/',map_name,'_AtlasLabel.dlabel.nii'))

#add in mwall mask values based on template
map$meta$cortex$medial_wall_mask <- mwall_template$meta$cortex$medial_wall_mask

#save out new cifti file
outfile <- paste0(dir,'/Final/',map_name,'_mwall')
write_cifti(map,outfile)