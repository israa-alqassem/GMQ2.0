
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""run_ca_dist.py: 
            run residuePair to get the distances between neigbouring residues
            """

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, GMQ"


# In[1]:

import sys,os
from auxiliary_function_and_variables import *


# In[5]:

pair_dist_dir = '/net/kihara/kang138/casp9-predictions/residuePair_1_2'
#target_list = get_pdb_file_paths().keys()
modellist = get_pdb_model_list()


# In[6]:

out_dir = '/net/kihara/ialqasse/Desktop/output/ca_dist/'

for model in modellist:
    
    print 'current model -->'+model
    
    os.chdir(out_dir)
    os.system('mkdir %s' % model)
    os.chdir(out_dir+model)
    
    target_pdbfile_dict = get_pdb_files(pdb_files_dir,model)
    
    for target in target_pdbfile_dict:
        pdbfile = target_pdbfile_dict[target]    
        os.system('%s/residuePair -f %s -t 6.0 -a > %s.ca.dist'                  %(pair_dist_dir, pdbfile, target))
    
    #os.system('mv %s.ca.dist ca_dist'%(target))


# In[ ]:



