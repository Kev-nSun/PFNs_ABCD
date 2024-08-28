import os
import sys
import time
import numpy as np
"""
syntax:  python submit_haufeAAB_1000_nulls.py mmddyy X (where X is A-J)

*****  this is only run network all null  *****

this is the only script you actually should call, it will submit a job
this job refers to the features, and then does the Haufe calculation. 
it assigns the current running dir to home_dir variable. It creates 
   tmp dir = <current dir>/tempfiles_matchedsamples_0405023 (based on already run ridge regression)
   result dir = <current dir>/Weight_Haufe_Nulls_<mmddyy>_All_cov_<X>
"""

dirdate = sys.argv[1]   # the command line parameter as mmddyy
jobID = sys.argv[2]     # the command line parameter as X (A-J)
homedir = os.getcwd()   # get the current running dir

tmpdir = '{0}/tempfiles_matchedsamples_0405023'.format(homedir)
sge_dir = '{0}/sge_{1}'.format(homedir,dirdate)  # dir for the processors std output and system error output files

os.makedirs(sge_dir,exist_ok=True)  # Create the output dir if it doesn't exist

networks = np.array(['all']) 
network = 'all'

# submit_ks.py has processed the preprocess_ks.py for all.  Here don't need to run again. If not, uncomment lines 32 & 33.
# following call preprocess_ks.py is moved from ks_proc_predict.py to avoid repeating for each psy_score run which causes the processor to die due to system resource issues.
#print('python {0}/preprocess_ks.py {1} {2}'.format(homedir,network,tmpdir))
#os.system('python {0}/preprocess_ks.py {1} {2}'.format(homedir,network,tmpdir))

GB = '400G'

for outcome in ['PRS_1','PRS_2']:
    print('qsub -l h_vmem={0},s_vmem={0} -N {1}{4} -R y -V -j y -b y -o ./sge_{6}/ -e ./sge_{6}/ python {2}/ks_proc_predict_haufeAAB_1000_nulls.py {2} {3} {4} {5} {6} {7}'.format(GB,outcome[1],homedir,tmpdir,network,outcome,dirdate,jobID))
    os.system('qsub -l h_vmem={0},s_vmem={0} -N {1}{4} -R y -V -j y -b y -o ./sge_{6}/ -e ./sge_{6}/ python {2}/ks_proc_predict_haufeAAB_1000_nulls.py {2} {3} {4} {5} {6} {7}'.format(GB,outcome[1],homedir,tmpdir,network,outcome,dirdate,jobID))
    #time.sleep(3)
