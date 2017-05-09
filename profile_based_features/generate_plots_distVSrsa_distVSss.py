
# coding: utf-8

# In[ ]:

"""generate_plots_distVSrsa_distVSss.py: 

    Generates plots of:
    
            distance vs. secondary structure agreement ratio
            distance vs.  RSA agreement ratio
            
    The distance is measured between the target and the model
    
"""

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, KiharaLab"


# In[1]:

import os
from pandas import *
import numpy as np
import itertools
import matplotlib.pyplot as plt
#from sklearn import svm, datasets
#from sklearn.metrics import confusion_matrix
from scipy.stats.stats import pearsonr   
from auxiliary_function_and_variables import *
get_ipython().magic(u'matplotlib inline')


# In[2]:

def file_exist(filename):
    if not os.path.exists(filename):
        return False

    if os.stat(filename).st_size == 0:
        return False
    
    return True


# In[3]:

# target folders -> model files (each file contains the distance between that model and the target)
calpha_dir = '/net/kihara/ialqasse/Desktop/CALPHA_result'
# model folders -> RSA_per_residue file for each target
rsa_results_dir = '/net/kihara/ialqasse/Desktop/output/RSA_results'
ss_results_dir = '/net/kihara/ialqasse/Desktop/output/SS_results'

RSA_figures_dir = '/net/kihara/ialqasse/Desktop/output/plots/RSA/'
SS_figures_dir = '/net/kihara/ialqasse/Desktop/output/plots/SS/'


# In[4]:

target_list = get_pdb_file_paths().keys()
model_list = get_pdb_model_list()


# In[5]:

dist_file_dict = {}
rsa_file_dict = {}
ss_file_dict = {}

for target in target_list:
    dist_file_dict[target] = []
    rsa_file_dict[target] = []
    ss_file_dict[target] = []


# In[6]:

for target in target_list:
    for model in model_list:
        
        if model.startswith('3D-JIGSAW'):
            continue
        
        distance_file = calpha_dir+'/'+target+'/'+model+'.lga'
        rsa_file = rsa_results_dir+'/'+model+'/'+target+'_RSA_per_residue.txt'
        ss_file = ss_results_dir+'/'+model+'/'+target+'_SS_per_residue.txt'
        
        # If the model file does not exist for any of ss/rsa/distance then we skip that model
        if(file_exist(distance_file) == False or file_exist(rsa_file) == False or file_exist(ss_file) == False):
            continue
        
        dist_file_dict[target].append(distance_file)
        rsa_file_dict[target].append(rsa_file)
        ss_file_dict[target].append(ss_file)


# In[7]:

def get_distance_list_per_residue(target):
    residue_distance_dict = {}
    for filename in dist_file_dict[target]:
        with open(filename) as f:
            next(f) # ignore header
             #lines = (line for line in lines if line) # Non-blank lines
            for line in f:
                if line.strip() != '':
                    line = line.split()
                    if line[0] == 'LGA':
                        model_RES = line[1]
                        real_RES = line[3]
                        id_RES = line[2].split('_')[0]
                        dist = line[5]
                        key = id_RES+'_'+model_RES

                        if key not in residue_distance_dict:
                            residue_distance_dict[key]=[]

                        residue_distance_dict[key].append(dist)
                        
    return residue_distance_dict


# In[8]:

def get_RSA_list_per_residue(target):
    residue_RSA_dict = {}
    for filename in rsa_file_dict[target]:
        if file_exist(filename) == False:
            continue
        with open(filename) as f:
            next(f)
            for line in f:
                line = line.split()
                res_id = line[0]
                res = line[1]
                RSA_agreement = line[2]
                key = res_id+'_'+res

                if key not in residue_RSA_dict:
                    residue_RSA_dict[key] = []

                residue_RSA_dict[key].append(RSA_agreement)
                
    return residue_RSA_dict


# In[9]:

def get_SS_list_per_residue(target):
    residue_SS_dict = {}
    for filename in ss_file_dict[target]:
        with open(filename) as f:
            if file_exist(filename) == False:
                continue
            header_line = next(f)
            for line in f:
                line = line.split()
                res_id = line[0]
                res = line[1]
                #real_class = line[2]
                #pred_class = line[3]
                SS_agreement = line[4]
                key = res_id+'_'+res

                if key not in residue_SS_dict:
                    residue_SS_dict[key] = []

                residue_SS_dict[key].append(SS_agreement)
                
    return residue_SS_dict


# In[29]:

target
print len(rsa_list)
print len(dist_list)


# In[22]:

#target = 'T0571'
for target in target_list: 
    residue_SS_dict = get_SS_list_per_residue(target)
    residue_RSA_dict = get_RSA_list_per_residue(target)
    residue_distance_dict = get_distance_list_per_residue(target)

    dist_list = []
    rsa_list = []
    ss_list = []

    for residue in residue_SS_dict.keys():
        if residue in residue_RSA_dict and residue in residue_distance_dict:
            dist_list.extend(residue_distance_dict[residue])
            rsa_list.extend(residue_RSA_dict[residue])
            ss_list.extend(residue_SS_dict[residue])


    ss_list = [float(i) for i in ss_list]  
    dist_list = [float(i) for i in dist_list]  
    rsa_list = [float(i) for i in rsa_list]  
    
    if len(dist_list) != 0 and len(rsa_list) != 0 and len(ss_list) != 0:
        
        print target
        print "Pearson correlation coefficient (RSA and distance) = ",
        print pearsonr(rsa_list,ss_list)
        print "\nPearson correlation coefficient (SS and distance) = ",
        print pearsonr(ss_list,dist_list)
    
        plt.figure()
        plt.xlabel('SS agreement ratio')
        plt.ylabel('Distance')
        plt.title(target)
        plt.plot(ss_list,dist_list,'.') 
        plt.savefig(SS_figures_dir+target+"2.png")

        plt.clf()

        plt.figure()
        plt.xlabel('RSA agreement ratio')
        plt.ylabel('Distance')
        plt.title(target)
        plt.plot(rsa_list,dist_list,'.') 
        plt.savefig(RSA_figures_dir+target+"2.png")

        plt.clf()

