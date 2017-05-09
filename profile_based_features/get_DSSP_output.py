
# coding: utf-8

# In[1]:

import sys,os
from auxiliary_function_and_variables import *


# In[4]:

def get_pdb_files(DIR,model_name,extension = '.mod'):

    # the immediate directories
    folder_list = [folder for folder in next(os.walk(DIR))[1]]
    target_pdbfile_dict = {}
    for folder in folder_list:
        if folder == 'modified_pdb':
            continue
        for f in os.listdir(DIR+'/'+folder):
            if f.endswith(extension) and f.startswith(model_name):
                target_pdbfile_dict[folder] = DIR+'/'+folder+'/'+f

    return target_pdbfile_dict


# In[7]:

#def main():

#os.system('mkdir dssp_target')
pdb_files_dir = '/net/kihara/ialqasse/Desktop/input_files/target_pdb'
dssp_dir = '/net/kihara/kang138/casp9-predictions/DSSP'
output_dir = '/net/kihara/ialqasse/Desktop/output/DSSP_output/'

target_list = get_pdb_file_paths().keys()
modellist = get_pdb_model_list()


for model in modellist:
    
    print 'current model -->'+model
    
    os.chdir(output_dir)
    os.system('mkdir %s' % model)

    target_pdbfile_dict = get_pdb_files(pdb_files_dir,model)
    os.chdir(output_dir+model)
    
    for target in target_pdbfile_dict:
        pdbfile = target_pdbfile_dict[target]    
        os.system('%s/dssp-2.0.4-linux-amd64 -i %s -o %s.dssp'                  %(dssp_dir, pdbfile, target))





# In[ ]:



