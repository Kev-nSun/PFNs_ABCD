import numpy as np
from scipy.io import savemat

# setup the data dir relative to the code running dir 
Datadir='Haufe_PRS_040724/Weight_Haufe_040724'

F1_A=np.load('{0}/PRS_1_A_Weight_haufetrans_all.npy'.format(Datadir))
mdic={"F1_A": F1_A}
savemat('{0}/PRS_1_A_Weight_haufetrans_all.mat'.format(Datadir),mdic)

F1_B=np.load('{0}/PRS_1_B_Weight_haufetrans_all.npy'.format(Datadir))
mdic={"F1_B": F1_B}
savemat('{0}/PRS_1_B_Weight_haufetrans_all.mat'.format(Datadir),mdic)

F2_A=np.load('{0}/PRS_2_A_Weight_haufetrans_all.npy'.format(Datadir))
mdic={"F2_A": F2_A}
savemat('{0}/PRS_2_A_Weight_haufetrans_all.mat'.format(Datadir),mdic)

F2_B=np.load('{0}/PRS_2_Weight_haufetrans_all.npy'.format(Datadir))
mdic={"F2_B": F2_B}
savemat('{0}/PRS_2_B_Weight_haufetrans_all.mat'.format(Datadir),mdic)

