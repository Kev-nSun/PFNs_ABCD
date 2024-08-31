#medial wall needs to be included with NA values, can use move_from_mwall(map)
#maps need to be in GIFTI format, can convert using C:\workbench\bin_windows64\ wb_command -cifti-convert -to-gifti-ext "input.nii" "output.gii"

from neuromaps.stats import compare_images
from neuromaps import nulls, datasets
import nibabel as nib
import numpy as np
import os

indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_6/BASELINE/weight_map_HaufeAAB/GIFTI"
GP_map_avg = nib.load("{0}/GPb_Haufe_AVG_FIRST.gii".format(indir))
GP_map_A = nib.load("{0}/GPb_Haufe_A_med.gii".format(indir))
GP_map_B = nib.load("{0}/GPb_Haufe_B_med.gii".format(indir))

indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/4_Brain_maps/FULL/weight_maps_HaufeAAB/GIFTI/"
F1_map_avg = nib.load("{0}F1_Haufe_AVG_FIRST.gii".format(indir))
F1_map_A = nib.load("{0}F1_Haufe_A_med.gii".format(indir))
F1_map_B = nib.load("{0}F1_Haufe_B_med.gii".format(indir))
F2_map_avg = nib.load("{0}F2_Haufe_AVG_FIRST.gii".format(indir))
F2_map_A = nib.load("{0}F2_Haufe_A_med.gii".format(indir))
F2_map_B = nib.load("{0}F2_Haufe_B_med.gii".format(indir))

outdir = 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/F1+F2_Scripts+files/Step_6/BASELINE-FULL/Spin_Results'
os.makedirs(outdir,exist_ok=True)  # Create the output dir if it doesn't exist

#spin nulls to compare between group A and B
GP_A_null = nulls.alexander_bloch(GP_map_A, atlas='fsLR', density='32k', n_perm=1000, seed=1234)
np.save('{0}/GP_A_null.npy'.format(outdir), GP_A_null)
F1_A_null = nulls.alexander_bloch(F1_map_A, atlas='fsLR', density='32k', n_perm=1000, seed=1234)
np.save('{0}/F1_A_null.npy'.format(outdir), F1_A_null)
F2_A_null = nulls.alexander_bloch(F2_map_A, atlas='fsLR', density='32k', n_perm=1000, seed=1234)
np.save('{0}/F2_A_null.npy'.format(outdir), GP_A_null)

#comparing between group A and B
corr_GP_A_B, p_GP_A_B = compare_images(GP_map_A, GP_map_B, nulls=GP_A_null)
corr_F1_A_B, p_F1_A_B = compare_images(F1_map_A, F1_map_B, nulls=F1_A_null)
corr_F2_A_B, p_F2_A_B = compare_images(F2_map_A, F2_map_B, nulls=F2_A_null)

spin_AB = np.array([[corr_GP_A_B, p_GP_A_B], [corr_F1_A_B, p_F1_A_B], [corr_F2_A_B, p_F2_A_B]])
np.savetxt('{0}/spin_AB_results.csv'.format(outdir), spin_AB, delimiter=",")

#spin nulls to compare between GP, F1, F2 avgs between A and B
GP_avg_null = nulls.alexander_bloch(GP_map_avg, atlas='fsLR', density='32k', n_perm=1000, seed=1234)
np.save('{0}/GP_avg_null.npy'.format(outdir), GP_avg_null)
F1_avg_null = nulls.alexander_bloch(F1_map_avg, atlas='fsLR', density='32k', n_perm=1000, seed=1234)
np.save('{0}/F1_avg_null.npy'.format(outdir), F1_avg_null)
F2_avg_null = nulls.alexander_bloch(F2_map_avg, atlas='fsLR', density='32k', n_perm=1000, seed=1234)
np.save('{0}/F2_avg_null.npy'.format(outdir), F2_avg_null)

#comparing between GP, F1, F2 avgs
corr_F1_F2, p_F1_F2 = compare_images(F1_map_avg, F2_map_avg, nulls=F1_avg_null)
corr_GP_F1, p_GP_F1 = compare_images(GP_map_avg, F1_map_avg, nulls=GP_avg_null)
corr_GP_F2, p_GP_F2 = compare_images(GP_map_avg, F2_map_avg, nulls=GP_avg_null)

spin_GP_F1_F2 = np.array([[corr_F1_F2, p_F1_F2], [corr_GP_F1, p_GP_F1], [corr_GP_F2, p_GP_F2]])
print(spin_GP_F1_F2)
np.savetxt('{0}/spin_GP_F1_F2_results.csv'.format(outdir), spin_GP_F1_F2, delimiter=",")

#comparing GP, F1, F2 avgs to SA axis
SAaxis = datasets.fetch_annotation(source='sydnor2021')
corr_GP_SA, p_GP_SA = compare_images(GP_map_avg, SAaxis, nulls=GP_avg_null)
corr_F1_SA, p_F1_SA = compare_images(F1_map_avg, SAaxis, nulls=F1_avg_null)
corr_F2_SA, p_F2_SA = compare_images(F2_map_avg, SAaxis, nulls=F2_avg_null)

spin_SA = np.array([[corr_GP_SA, p_GP_SA], [corr_F1_SA, p_F1_SA], [corr_F2_SA, p_F2_SA]])
print(spin_SA)
np.savetxt('{0}/spin_SA_results.csv'.format(outdir), spin_SA, delimiter=",")


indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_6/weight_map_by_vertex_HaufeAAB_All_cov_nii/GIFTI"
GP1_map_avg = nib.load("{0}/PB1_Haufe_AVG_FIRST.gii".format(indir))

#comparing between GP baseline and y1
corr_GP_long, p_GP_long = compare_images(GP_map_avg, GP1_map_avg, nulls=GP_avg_null)
spin_GP_long = np.array([corr_GP_long, p_GP_long])
print(spin_GP_long)
np.savetxt('{0}/spin_GP_long_results.csv'.format(outdir), spin_GP_long, delimiter=",")