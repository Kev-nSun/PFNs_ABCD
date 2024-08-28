import sys
import os

homedir = sys.argv[1]  #'/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples'
tmpdir = sys.argv[2]   #'/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples/tempfiles_matchedsamples_mmddyy'
network = sys.argv[3]  # all
phenotype_name = sys.argv[4]   #General_PB1
dirdate = sys.argv[5]   # the date as mmddyy
jobID = sys.argv[6]

# the following call is moved to submit_ks.py to avoid repeating for each psy_score run which causes the processor to die due to system resource issues.
#print('python {0}/preprocess_ks.py {1} {2}'.format(homedir,network,tmpdir))
#os.system('python {0}/preprocess_ks.py {1} {2}'.format(homedir,network,tmpdir))

output_dir = homedir + '/Weight_Haufe_Nulls_{0}_All_cov_{1}'.format(dirdate,jobID)
feature_path = '{0}/features_{1}.npy'.format(tmpdir,network)
phenotype_path = '{0}/phenotypes_target.csv'.format(tmpdir)
pred_path_dir = homedir + '/results_matchedsamples_null_1000_020924' + jobID

"""
this does the actual Haufe calculation
"""
print('python {0}/Weight_haufeAAB_1000_nulls.py {1} {2} {3} {4} {5} rel_family_id'.format(homedir,output_dir,feature_path,phenotype_path,phenotype_name,pred_path_dir))
os.system('python {0}/Weight_haufeAAB_1000_nulls.py {1} {2} {3} {4} {5} rel_family_id'.format(homedir,output_dir,feature_path,phenotype_path,phenotype_name,pred_path_dir))
