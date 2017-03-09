
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""run_naccess_and_ACCPro.py: 
        generates output files for actual and predicted surface areas
            a. runs naccess and stores the output files in folder (whose name is the same as target name/ID) .
            b. runs ACCPro and keeps the output in ACCPro_output directory  """


# In[17]:

import os
from auxiliary_function_and_variables import *


# In[18]:

target_files_dict = get_fasta_pdb_dict()
NACCESS_output = '/net/kihara/ialqasse/Desktop/output/NACCESS_output'
SCRATCH_DIR = '/net/kihara/ialqasse/SCRATCH-1D_1.1/bin'
ACCPro_output = '/net/kihara/ialqasse/Desktop/output/ACCPro_output'


# In[19]:

for target in target_files_dict.keys():
    # Actual surface area
    os.chdir(NACCESS_output)
    os.system('mkdir %s' %(target))
    os.chdir(NACCESS_output+'/'+target)
    os.system('naccess %s' %(target_files_dict[target][1]))
    
    # Predicted surface area
    os.chdir(ACCPro_output)
    os.system('%s/run_SCRATCH-1D_predictors.sh %s %s 4 ' %(SCRATCH_DIR, target_files_dict[target][0], target))
    

