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
This script will run the Haufe transform (AAB) on the weights from a previously run ridge regression of all PFNs.

"""


output_dir = sys.argv[1] #'./Weight_Haufe_Nulls_mmddyy_X'
os.makedirs(output_dir,exist_ok=True)
feature_path = sys.argv[2]  #'../Step_2/120222_tempfiles/tempfiles_matchedsamples_031324/features_all.npy'
phenotype_path = sys.argv[3] #'../Step_2/120222_tempfiles/tempfiles_matchedsamples_031324/phenotypes_target.csv'
phenotype_name = sys.argv[4] #'General_PB1'
fold_group = 'matched_group' #sys.argv[6] # in this instance, this is the ABCD train and test split, A for one group, B for the other
pred_path_dir = sys.argv[5] # '/results_matchedsamples_null_1000_031424X'


print(output_dir,end='\n')
print(feature_path,end='\n')
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

var_A = np.var(targets[A]) # variance of actual p-factor in group A
var_B = np.var(targets[B]) # variance of actual p-factor in group B
print (var_A)
print (var_B)

features = np.load(feature_path).astype(np.float16)
print('load ',feature_path,'\n')
#print("features type: ",type(features),end='\n')
#print(np.shape(features),end='\n')

"""
this is adapted from pennlinckit.utils.predict
"""

# np.save('{0}/{1}_targets.npy'.format(output_dir,phenotype_name),targets)
#assert targets.shape[0] == features.shape[0]

#accuracy = []
#prediction = np.zeros((targets.shape))
feat_A = np.matrix(features[A]).tolist()
print("A: x type: ",type(feat_A),end='\n')
#print(np.shape(x),end='\n')

feat_B = np.matrix(features[B]).tolist()
print("B: x type: ",type(feat_B),end='\n')
#print(np.shape(x),end='\n')

for fold in range(100):
    #------------------------ for A -------------------------------
    print('loop',fold)

    pred = np.load('{0}/all_network/{1}_prediction_testA_{2}.npy'.format(pred_path_dir,phenotype_name,fold))
    print('load {0}/all_network/{1}_prediction_testA_{2}.npy'.format(pred_path_dir,phenotype_name,fold))
    #print("pred: ",type(pred),end='\n')
    #print(np.shape(pred),end='\n')

    haufetrans = []

    # for subjects in group A, calculate covariance between ridge-predicted p-factor and loading value (one vertex at a time)
    for vert in range(1010004):
        covmat = np.cov(pred, [sub[vert] for sub in feat_A]) # returns column corresponding to 'vert' in 'feat_A' (loading values)
        cov = covmat[0,1]
        if math.isnan(cov): # if cov is Na, replace with 0
            cov = 0.0
        cov_norm = cov/var_A # normalize by variance of actual p-factor in group A
        haufetrans.append(cov_norm)

    #print("haufetrans: ",type(haufetrans),np.shape(haufetrans),end='\n')
    #print(np.shape(haufetrans),end='\n')

    print('haufetrans save to {0}/{1}_A_Weight_haufetrans_all_{2}.npy'.format(output_dir,phenotype_name,fold))
    np.save('{0}/{1}_A_Weight_haufetrans_all_{2}.npy'.format(output_dir,phenotype_name,fold),haufetrans)

    #------------------------ for B -------------------------------
    pred = np.load('{0}/all_network/{1}_prediction_testB_{2}.npy'.format(pred_path_dir,phenotype_name,fold))
    print('load {0}/all_network/{1}_prediction_testB_{2}.npy'.format(pred_path_dir,phenotype_name,fold))
    #print("pred: ",type(pred),end='\n')
    #print(np.shape(pred),end='\n')

    haufetrans = []

    # for subjects in group B, calculate covariance between ridge-predicted p-factor and loading value (one vertex at a time)
    for vert in range(1010004):
        covmat = np.cov(pred, [sub[vert] for sub in feat_B]) # sub[vert] returns column corresponding to 'vert' in 'feat_B' (loading values)
        cov = covmat[0,1]
        if math.isnan(cov): # if cov is Na, replace with 0
            cov = 0.0
        cov_norm = cov/var_B # normalize by variance of actual p-factor in group B
        haufetrans.append(cov_norm)

    #print("haufetrans: ",type(haufetrans),np.shape(haufetrans),end='\n')
    #print(np.shape(haufetrans),end='\n')
    print('haufetrans save to {0}/{1}_B_Weight_haufetrans_all_{2}.npy'.format(output_dir,phenotype_name,fold))
    np.save('{0}/{1}_B_Weight_haufetrans_all_{2}.npy'.format(output_dir,phenotype_name,fold),haufetrans)
