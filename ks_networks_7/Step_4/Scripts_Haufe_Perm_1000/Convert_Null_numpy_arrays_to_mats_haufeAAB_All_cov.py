import numpy as np
import sys
from scipy.io import savemat

# call script: 'python Convert_Null_numpy_arrays_to_mats_haufeAAB_All_cov.py X', where X is a letter A-J
# setup the data dir relative to the code running dir
Dir=sys.argv[1]
Datadir='../Haufe_Nulls_031824/Weight_Haufe_Nulls_031824_{0}'.format(Dir)

#Cycle through and convert 200 npy files (100 A and 100 B) in each dir
for fold in range(100):
	PB_test=np.load('{0}/General_PB1_A_Weight_haufetrans_all_{1}.npy'.format(Datadir,fold))
	mdic={"PB_test": PB_test}
	savemat('{0}/General_PB1_A_Weight_haufetrans_all_{1}.mat'.format(Datadir,fold),mdic)

	PB_test=np.load('{0}/General_PB1_B_Weight_haufetrans_all_{1}.npy'.format(Datadir,fold))
	mdic={"PB_test": PB_test}
	savemat('{0}/General_PB1_B_Weight_haufetrans_all_{1}.mat'.format(Datadir,fold),mdic)
