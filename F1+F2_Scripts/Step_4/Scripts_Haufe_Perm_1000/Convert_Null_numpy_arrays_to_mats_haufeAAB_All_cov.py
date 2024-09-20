import numpy as np
import sys
from scipy.io import savemat

# call script: 'python Convert_Null_numpy_arrays_to_mats_haufeAAB_All_cov.py X', where X is a letter A-J
# setup the data dir relative to the code running dir
Dir=sys.argv[1]
Datadir='Weight_Haufe_Nulls_041024/Weight_Haufe_Nulls_041024_All_cov_{0}'.format(Dir)

#Cycle through anc convert 200 npy files (100 A and 100 B) in each dir

for PRS in ['PRS_1','PRS_2']:
	for fold in range(100):
		Weights=np.load('{0}/{1}_A_Weight_haufetrans_all_{2}.npy'.format(Datadir,PRS,fold))
		mdic={"Weights": Weights}
		savemat('{0}/{1}_A_Weight_haufetrans_all_{2}.mat'.format(Datadir,PRS,fold),mdic)

		Weights=np.load('{0}/{1}_B_Weight_haufetrans_all_{2}.npy'.format(Datadir,PRS,fold))
		mdic={"Weights": Weights}
		savemat('{0}/{1}_B_Weight_haufetrans_all_{2}.mat'.format(Datadir,PRS,fold),mdic)
