
# coding: utf-8

# In[1]:

import os


# In[7]:

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


# In[6]:

def get_pdb_model_list(DIR):
    extension = '.mod' # same as '.pdb'
    # the immediate directories
    folder_list = [folder for folder in next(os.walk(DIR))[1]]
    for folder in folder_list:
        if folder == 'modified_pdb':
            continue

        modellist = os.listdir(DIR+'/'+folder)
                
        modellist = [modellist[i].split('.')[0] for i in range(0,len(modellist))]
                
        return modellist


# In[1]:

STRIDE_DIR = '/net/kihara/ialqasse/STRIDE'
pdb_files_dir = '/net/kihara/ialqasse/Desktop/input_files/target_pdb'
output_DIR = '/net/kihara/ialqasse/Desktop/output/pdb_models_output/'

modellist = get_pdb_model_list(pdb_files_dir)

for model in modellist:
    
    print 'current model -->'+model
    
    
    os.chdir(output_DIR)
    os.system('mkdir %s' % model)

    target_pdbfile_dict = get_pdb_files(pdb_files_dir,model)

    for target in target_pdbfile_dict:
        pdbfile = target_pdbfile_dict[target]
        
        os.chdir(output_DIR+model)
        os.system('mkdir %s' % target)
        os.chdir(output_DIR+model+'/'+target)
        
        os.system('%s/stride %s > %s' %(STRIDE_DIR,pdbfile,target))
        os.system('naccess %s' %pdbfile)
        


# In[ ]:



