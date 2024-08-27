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

Command line:

python Weight_haufeAAB_All_cov.py 
"""


output_dir = 'Weight_Haufe_031324_All_cov'
os.makedirs(output_dir,exist_ok=True)
feature_path_dir = 'tempfiles_matchedsamples_031324/'
phenotype_path = 'tempfiles_matchedsamples_031324/phenotypes_target.csv'
phenotype_name = 'General_PB1'
fold_group = 'matched_group' # in this instance, this is the ABCD train and test split, A for one group, B for the other
pred_path_dir = 'results_matchedsamples_031324/'


print(output_dir,end='\n')
print(feature_path_dir,end='\n')
print(phenotype_path,end='\n')
print(phenotype_name,end='\n')
print(pred_path_dir,end='\n')


targets = pd.read_csv(phenotype_path)[phenotype_name].values.astype(np.float16)
print('load targets - {0}'.format(phenotype_path))
#print("targets type: ",type(targets),end='\n')
#print(np.shape(targets),end='\n')
fold_group = pd.read_csv(phenotype_path)[fold_group].values.astype(np.float16)
#print('load fold_group -- {0}'.format(phenotype_path))
#print(fold_group,end='\n')
#print("fold_group type: ",type(fold_group),end='\n')
#print(np.shape(fold_group),end='\n')

# Split up the data into train/test based on the fold group
A = np.argwhere(fold_group==1) # <- save an index of everywhere the fold group is 1
B = np.argwhere(fold_group==2) # <- save an index of everywhere the fold group is 2
	
	
print("A type: ",type(A),np.shape(A),end='\n')
print("B type: ",type(B),np.shape(B),end='\n')

var_A = np.var(targets[A])
var_B = np.var(targets[B])
print (var_A)
print (var_B)

n = 'all'
feature_path = '{0}features_{1}.npy'.format(feature_path_dir,n)
features = np.load(feature_path).astype(np.float16)
print('load ',feature_path,'\n')
#print("features type: ",type(features),end='\n')
#print(np.shape(features),end='\n')


#------------------------ for A -------------------------------
x = np.matrix(features[A]).tolist()
print("A: x type: ",type(x),end='\n')
#print(np.shape(x),end='\n')

#extract feature importance
pred = np.load('{0}{1}_network/{2}_prediction_testA.npy'.format(pred_path_dir,n,phenotype_name))
print('load {0}{1}_network/{2}_prediction_testA.npy'.format(pred_path_dir,n,phenotype_name))
#print("pred: ",type(pred),end='\n')
#print(np.shape(pred),end='\n')

haufetrans = []

for i in range(1010004):
    covmat = np.cov(pred, [sub[i] for sub in x])
    cov = covmat[0,1]
    if math.isnan(cov):
        cov = 0.0
    cov_norm = cov/var_A
    haufetrans.append(cov_norm)

#print("haufetrans: ",type(haufetrans),np.shape(haufetrans),end='\n')
#print(np.shape(haufetrans),end='\n')

print('haufetrans save to {0}/{1}_A_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n))
np.save('{0}/{1}_A_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n),haufetrans)

#------------------------ for B -------------------------------
x = np.matrix(features[B]).tolist()
print("B x type: ",type(x),end='\n')
#print(np.shape(x),end='\n')

#extract feature importance
pred = np.load('{0}{1}_network/{2}_prediction_testB.npy'.format(pred_path_dir,n,phenotype_name))
print('load {0}{1}_network/{2}_prediction_testB.npy'.format(pred_path_dir,n,phenotype_name))
#print("pred: ",type(pred),end='\n')
#print(np.shape(pred),end='\n')

haufetrans = []

for i in range(1010004):
    covmat = np.cov(pred, [sub[i] for sub in x])
    cov = covmat[0,1]
    if math.isnan(cov):
        cov = 0.0
    cov_norm = cov/var_B
    haufetrans.append(cov_norm)

#print("haufetrans: ",type(haufetrans),np.shape(haufetrans),end='\n')
#print(np.shape(haufetrans),end='\n')
print('haufetrans save to {0}/{1}_B_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n) )
np.save('{0}/{1}_B_Weight_haufetrans_{2}.npy'.format(output_dir,phenotype_name,n),haufetrans)
  
  