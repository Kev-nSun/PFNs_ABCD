#maps need to be in GIFTI format, can convert using C:\workbench\bin_windows64\wb_command -cifti-convert -to-gifti-ext "input.nii" "output.gii"
#used bat script "convert.bat" to cycle through 1000 ciftis to convert to giftis

from neuromaps.stats import compare_images
from neuromaps import nulls, datasets
import nibabel as nib
import numpy as np
import os
from scipy.io import savemat

GP_F1 = []
GP_F2 = []
GP_SA = []
GP_AB = []

F1_F2 = []
F1_SA = []
F1_AB = []

F2_SA = []
F2_AB = []

indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_6/BASELINE/weight_map_HaufeAAB/GIFTI/" #dir for GP actual
GP_map_avg = nib.load("{0}GPb_Haufe_AVG_FIRST.gii".format(indir))
GP_map_A = nib.load("{0}GPb_Haufe_A_med.gii".format(indir))
GP_map_B = nib.load("{0}GPb_Haufe_B_med.gii".format(indir))

indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/4_Brain_maps/FULL/weight_maps_HaufeAAB/GIFTI/" #dir for PRS actual
F1_map_avg = nib.load("{0}F1_Haufe_AVG_FIRST.gii".format(indir))
F1_map_A = nib.load("{0}F1_Haufe_A_med.gii".format(indir))
F1_map_B = nib.load("{0}F1_Haufe_B_med.gii".format(indir))
F2_map_avg = nib.load("{0}F2_Haufe_AVG_FIRST.gii".format(indir))
F2_map_A = nib.load("{0}F2_Haufe_A_med.gii".format(indir))
F2_map_B = nib.load("{0}F2_Haufe_B_med.gii".format(indir))

SAaxis = datasets.fetch_annotation(source='sydnor2021')

for fold in range(1000): #loop through null maps and correlate with actual maps
	num = fold + 1 #file numbers start from 1--> 1000 (loop starts at 0)	

	#correlating null GP avg to F1, F2 avgs and SA
	indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_6/BASELINE/weight_map_HaufeAAB_nulls/GIFTI/" #dir for GP nulls
	Null_GP_map_avg = nib.load("{0}/General_PB1_Avg_Weight_haufetrans_all_med_wall_{1}.gii".format(indir,num)) #access GP avg null gifti
	corr = compare_images(Null_GP_map_avg, F1_map_avg, metric='pearsonr')
	GP_F1.append(corr) 
	corr = compare_images(Null_GP_map_avg, F2_map_avg, metric='pearsonr')
	GP_F2.append(corr)
	corr = compare_images(Null_GP_map_avg, SAaxis, metric='pearsonr')
	GP_SA.append(corr)
	#correlating actual GP A to null B
	Null_GP_map_B = nib.load("{0}/General_PB1_B_Weight_haufetrans_all_med_wall_{1}.gii".format(indir,num)) #access GP B null gifti
	corr = compare_images(GP_map_A, Null_GP_map_B, metric='pearsonr')
	GP_AB.append(corr) 

	#correlating null F1 avg to F2 avg and SA
	indir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/F1+F2_Scripts+files/Step_5/FULL/Null_weight_maps/GIFTI/" #set indir to F1, F2 nulls
	Null_F1_map_avg = nib.load("{0}/PRS_1_Avg_Weight_haufetrans_all_med_wall_{1}.gii".format(indir,num)) #access F1 avg null gifti
	corr = compare_images(Null_F1_map_avg, F2_map_avg, metric='pearsonr')
	F1_F2.append(corr)
	corr = compare_images(Null_F1_map_avg, SAaxis, metric='pearsonr')
	F1_SA.append(corr)
	#correlating actual F1 A to null B
	Null_F1_map_B = nib.load("{0}/PRS_1_B_Weight_haufetrans_all_med_wall_{1}.gii".format(indir,num)) #access F1 B null gifti
	corr = compare_images(F1_map_A, Null_F1_map_B, metric='pearsonr')
	F1_AB.append(corr)

	#correlating null F2 avg to SA
	Null_F2_map_avg = nib.load("{0}/PRS_2_Avg_Weight_haufetrans_all_med_wall_{1}.gii".format(indir,num)) #access F2 avg null gifti
	corr = compare_images(Null_F2_map_avg, SAaxis, metric='pearsonr')
	F2_SA.append(corr)
	#correlating actual F2 A to null B
	Null_F2_map_B = nib.load("{0}/PRS_2_B_Weight_haufetrans_all_med_wall_{1}.gii".format(indir,num)) #access F2 B null gifti
	corr = compare_images(F2_map_A, Null_F2_map_B, metric='pearsonr')
	F2_AB.append(corr)



outdir = 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/F1+F2_Scripts+files/Step_6/BASELINE-FULL/Spin_Results'
os.makedirs(outdir,exist_ok=True)  # Create the output dir if it doesn't exist

mdic = {"GP_F1": GP_F1, "GP_F2": GP_F2, "GP_SA": GP_SA}
savemat("{0}/Corr_GP_nulls_avg.mat".format(outdir), mdic)
mdic = {"GP_AB": GP_AB}
savemat("{0}/Corr_GP_nulls_AB.mat".format(outdir), mdic)

mdic = {"F1_F2": F1_F2, "F1_SA": F1_SA}
savemat("{0}/Corr_F1_nulls_avg.mat".format(outdir), mdic)
mdic = {"F1_AB": F1_AB}
savemat("{0}/Corr_F1_nulls_AB.mat".format(outdir), mdic)

mdic = {"F2_SA": F2_SA}
savemat("{0}/Corr_F2_nulls_avg_SA.mat".format(outdir), mdic)
mdic = {"F2_AB": F2_AB}
savemat("{0}/Corr_F2_nulls_AB.mat".format(outdir), mdic)

