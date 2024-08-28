import glob
import h5py
import numpy as np
import os
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
import sys
from scipy.stats import pearsonr
import math

"""
This script will run the Haufe transform (AAB) on the weights from a previously run ridge regression of all PFNs. To run it, you need:

output_dir: where you want your results saved, does not have to exist, we will make it 
feature_path: path to the actual features, should be a subject x features npy array
phenotype_path: path to a csv with the to-be-predicted values
phenotype_name: str of name of column in phenotype_path you want
control_path: path to csv, regress these features from the features within each fold
fold_group: in this instance, this is the ABCD train and test split, A for one group, B for the other

"""


output_dir = sys.argv[1] #'Weight_Haufe_mmddyy'
os.makedirs(output_dir,exist_ok=True)
feature_path = sys.argv[2] #'tempfiles_matchedsamples_040524/'
phenotype_path = sys.argv[3] #'tempfiles_matchedsamples_040524/phenotypes_target.csv'
phenotype_name = sys.argv[4] #'PRS_1' or 'PRS_2'
fold_group = 'matched_group' # in this instance, this is the ABCD train and test split, A for one group, B for the other
pred_path_dir = sys.argv[5] #'results_matchedsamples_040524/'


# load in phenotypes (PRS_1 or PRS_2)
targets = pd.read_csv(phenotype_path)[phenotype_name].values.astype(np.float16)
print('load targets - {0}'.format(phenotype_path))

#load in fold group
fold_group = pd.read_csv(phenotype_path)[fold_group].values.astype(np.float16)

# Split up the data into train/test based on the fold group
A = np.argwhere(fold_group==1) # <- save an index of everywhere the fold group is 1
B = np.argwhere(fold_group==2) # <- save an index of everywhere the fold group is 2
	
print("A type: ",type(A),np.shape(A),end='\n')
print("B type: ",type(B),np.shape(B),end='\n')

#calculate variance of phenotype in each fold group
var_A = np.var(targets[A])
var_B = np.var(targets[B])
print (var_A)
print (var_B)

n = 'all'
features = np.load(feature_path).astype(np.float16)
#print("features type: ",type(features),end='\n')
#print(np.shape(features),end='\n')


#------------------------ for A -------------------------------
x = np.matrix(features[A]).tolist()
print("A: x type: ",type(x),end='\n')
#print(np.shape(x),end='\n')

#extract feature importance
pred = np.load('{0}{1}_network/{2}_prediction_testA.npy'.format(pred_path_dir,n,phenotype_name))
print('load {0}{1}_network/{2}_prediction_testA.npy'.format(pred_path_dir,n,phenotype_name))

haufetrans = [] # empty matrix

for i in range(1010004):
    covmat = np.cov(pred, [sub[i] for sub in x])
    cov = covmat[0,1]
    if math.isnan(cov):
        cov = 0.0
    cov_norm = cov/var_A
    haufetrans.append(cov_norm)


print('haufetrans save to {0}/{1}_A_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n))
np.save('{0}/{1}_A_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n),haufetrans)

#------------------------ for B -------------------------------
x = np.matrix(features[B]).tolist()
print("B x type: ",type(x),end='\n')
#print(np.shape(x),end='\n')

#extract feature importance
pred = np.load('{0}{1}_network/{2}_prediction_testB.npy'.format(pred_path_dir,n,phenotype_name))
print('load {0}{1}_network/{2}_prediction_testB.npy'.format(pred_path_dir,n,phenotype_name))

haufetrans = [] # empty matrix

for i in range(1010004):
    covmat = np.cov(pred, [sub[i] for sub in x])
    cov = covmat[0,1]
    if math.isnan(cov):
        cov = 0.0
    cov_norm = cov/var_B
    haufetrans.append(cov_norm)

print('haufetrans save to {0}/{1}_B_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n) )
np.save('{0}/{1}_B_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n),haufetrans)
  
  