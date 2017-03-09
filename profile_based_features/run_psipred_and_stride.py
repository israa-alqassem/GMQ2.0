
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""run_psipred_and_stride.py: 
        generates output files for actual and predicted secondary structure
            a. runs STRIDE and stores the output files in STRIDE_output directory .
            b. runs ACCPro and keeps the output in ACCPro_output directory  """


# In[1]:

import os
from auxiliary_function_and_variables import *


# In[4]:

target_file_dict = get_fasta_pdb_dict()
PSIPRED_DIR = '/net/kihara/ialqasse/PSIPRED/psipred'
STRIDE_DIR = '/net/kihara/ialqasse/STRIDE'
psipred_output_dir = '/net/kihara/ialqasse/Desktop/output/PSIPRED_output'
stride_output_dir = '/net/kihara/ialqasse/Desktop/output/STRIDE_output'


# In[6]:

for target in target_files_dict:
    fastafile = target_file_dict[target][0]
    pdbfile = target_file_dict[target][1]
    
    os.chdir(psipred_output_dir)
    os.system('%s/runpsipred %s' %(PSIPRED_DIR,fastafile))
    
    os.chdir(stride_output_dir)
    os.system('%s/stride %s > %s' %(STRIDE_DIR,pdbfile,target))

