import glob
import h5py
import numpy as np
import os
import pandas as pd
import sys
import scipy.stats
from scipy.io import savemat

"""
This script will calculage the mean absolute deviation (MAD) of each loading value 1010004 (17 x 59,412) across 7459 participants in ABCD baseline data, using 'features.npy' file produced during a prior ridge
regression code run, housed in  "tempfiles_matchedsamples_031324"

Will produce a mat file of 1010004 MAD values corresponding to cortical vertices

output_dir: where you want your results saved, does not have to exist, we will make it 
feature_path: path to the actual features, should be a subject x features npy array

Command line:

python PFN_loading_MAD_var.py PFN_loading_MAD_060324/ tempfiles_matchedsamples_031324/

"""

output_dir = sys.argv[1] #'PFN_loading_MAD_060324/'
os.makedirs(output_dir,exist_ok=True)
feature_path_dir = sys.argv[2]  #'tempfiles_matchedsamples_031324/'

feature_path = '{0}features_all.npy'.format(feature_path_dir)
features = np.load(feature_path).astype(np.float16)
print(features.shape)

MAD = []
for loading in range(1010004): #loop through each loading and derive MAD for each
    MAD_tracker = scipy.stats.median_abs_deviation(features[:,loading])
    MAD.append(MAD_tracker)

mdic={"MAD": MAD}
savemat('{0}PFN_loading_MAD.mat'.format(output_dir),mdic)