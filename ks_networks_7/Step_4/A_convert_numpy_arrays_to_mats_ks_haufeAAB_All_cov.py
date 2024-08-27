import numpy as np
from scipy.io import savemat

# setup the data dir relative to the code running dir 
Datadir='Weight_Haufe_031324'

PB_test=np.load('{0}/General_PB1_A_Weight_haufetrans_all.npy'.format(Datadir))
mdic={"PB_test": PB_test}
savemat('{0}/General_PB1_A_Weight_haufetrans_all.mat'.format(Datadir),mdic)

PB_test=np.load('{0}/General_PB1_B_Weight_haufetrans_all.npy'.format(Datadir))
mdic={"PB_test": PB_test}
savemat('{0}/General_PB1_B_Weight_haufetrans_all.mat'.format(Datadir),mdic)

