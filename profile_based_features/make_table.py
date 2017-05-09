
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""make_table.py: 
        This code generates training and test files used later on as input to prepare_gmq_training_file, prepare_gmq_test_file
        Each row of the output files has the following format:
        target~model residue_id distance cons ss_agreement rsa_agreement overall_rsa_agreement overall_ss_agreement
           
           
"""

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, GMQ"


# In[10]:

import sys,os
import numpy
from auxiliary_function_and_variables import *
import random


# In[11]:

ss_out_dir = '/net/kihara/ialqasse/Desktop/output/SS_results/'
rsa_out_dir = '/net/kihara/ialqasse/Desktop/output/RSA_results/'
psiblast_out_dir ='/net/kihara/ialqasse/Desktop/output/PSIBlast_output/'
calpha_dir = '/net/kihara/ialqasse/Desktop/CALPHA_result'


# In[12]:

def get_pdb_files(model_name):
    DIR = '/net/kihara/ialqasse/Desktop/input_files/target_pdb'
    extension = '.mod'
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


# In[13]:

def get_distance_list_per_residue(model,target):
    filename = calpha_dir+'/'+target+'/'+model+'.lga'
    residue_distance_dict = {}
    
    mean = 0.
    
    if not file_exist(filename):
        return residue_distance_dict,mean
    
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
                    dist = float(line[5])
                    key = id_RES+'_'+model_RES
                    residue_distance_dict[key] = dist
                    mean += dist
                        
    mean /= float(len(residue_distance_dict))
    return residue_distance_dict, mean


# In[14]:

def get_cons_per_residue(target):
    filename = psiblast_out_dir+target+'_PSSM'
    cons_res_dict = {}
    mean = 0.
    
    if not file_exist(filename):
        return cons_res_dict,mean
    

    with open(filename) as f:
        next(f) # ignore header lines
        next(f)
        next(f)
        for line in f:
            if line.strip() != '':
                line = line.split()
                res_id = line[0]

                # we read data of all residues
                if not res_id.isdigit():
                    break

                res = line[1]
                cons = float(line[42])
                #print res_id,
                #print res,
                #print cons
                cons_res_dict[res_id+'_'+res] = cons
                mean += cons
                
    mean /= float(len(cons_res_dict))
    return cons_res_dict,mean


# In[15]:

def get_ss_agreement_per_residue(model,target):
    
    ss_agreement_res_dict = {}
    mean = 0.
    
    filename = ss_out_dir+'/'+model+'/'+target+'_SS_per_residue.txt'
    
    if not file_exist(filename):
        return ss_agreement_res_dict,mean

    with open(filename) as f:
        next(f) # ignore header line
        for line in f:
            line = line.split()
            res_id = line[0]
            res = line[1]
            ss_agreement = float(line[4])
            ss_agreement_res_dict[res_id+'_'+res] = ss_agreement
            mean += ss_agreement
    mean /= float(len(ss_agreement_res_dict))  
    return ss_agreement_res_dict, mean


# In[16]:

def get_rsa_agreement_per_residue(model,target):
    
    rsa_agreement_res_dict = {}
    mean = 0.
    
    filename = rsa_out_dir+'/'+model+'/'+target+'_RSA_per_residue.txt'
    
    if not file_exist(filename):
        return rsa_agreement_res_dict,mean

    with open(filename) as f:
        next(f) # ignore header line
        for line in f:
            line = line.split()
            res_id = line[0]
            res = line[1]
            rsa_agreement = float(line[2])
            rsa_agreement_res_dict[res_id+'_'+res] = rsa_agreement
            mean += rsa_agreement
            
    mean /= float(len(rsa_agreement_res_dict))
    return rsa_agreement_res_dict, mean


# In[17]:

def get_rsa_agreement_per_target(model):
    
    target_rsa_dict = {}
    
    filename = rsa_out_dir+'/'+model+'/'+'overall_RSA.txt'
    
    if not file_exist(filename):
        return target_rsa_dict

    with open(filename) as f:
        next(f) # ignore header line
        for line in f:
            line = line.split()
            target = line[0]
            target_rsa_dict[target] = float(line[1])
            
    return target_rsa_dict


# In[18]:

def get_ss_agreement_per_target(model):
    
    target_ss_dict = {}
    
    filename = ss_out_dir+'/'+model+'/'+'overall_SS.txt'
    
    if not file_exist(filename):
        return target_rsa_dict
    

    with open(filename) as f:
        next(f) # ignore header line
        for line in f:
            line = line.split()
            target = line[0]
            target_ss_dict[target] = float(line[1])
            
    return target_ss_dict


# In[19]:

target_list = get_pdb_file_paths().keys()
modellist = get_pdb_model_list()

feature_list = []
train_set = []
test_set = []

train_flag = True

count = 0

for model in modellist:
    
    if model.startswith('3D-JIGSAW_'):
        continue
    
    target_rsa_dict = get_rsa_agreement_per_target(model)
    target_ss_dict = get_ss_agreement_per_target(model)
    target_list = get_pdb_files(model).keys()
    
    overall_mean_rsa = numpy.mean(target_rsa_dict.values()) ## overall mean for the whole model
    overall_mean_ss = numpy.mean(target_ss_dict.values()) 

    
    for target in target_list:
        #target = 'T0519'
        # Per target Split data into training and testing 8:2
        random_num = random.random()
        if random_num <= 0.20:
            train_flag = False
        else:
            train_flag = True


        res_cons_dict, mean_cons = get_cons_per_residue(target)
        ss_agreement_res_dict, mean_ss = get_ss_agreement_per_residue(model,target)
        rsa_agreement_res_dict, mean_rsa = get_rsa_agreement_per_residue(model,target)
        distance_res_dict, mean_dist = get_distance_list_per_residue(model,target)

        
        len1 = len(res_cons_dict)
        len2 = len(ss_agreement_res_dict)
        len3 = len(rsa_agreement_res_dict)
        
        if len1 != len2 or len1 != len3:
            continue
        else:
            count += 1
        
        # if we don't have the distance values, then skip it
        if len(distance_res_dict) == 0:
            continue

        if target in target_rsa_dict:
            overall_rsa = target_rsa_dict[target]
        else:
            overall_rsa = overall_mean_rsa
            #continue

        if target in target_ss_dict:
            overall_ss = target_ss_dict[target]
        else:
            overall_ss = overall_mean_ss
            #continue

        for key in ss_agreement_res_dict.keys():
            res_id = int(key.split('_')[0])
            #cons = res_cons_dict[key]

            if key in ss_agreement_res_dict:
                ss_agreement = ss_agreement_res_dict[key]
            else:
                ss_agreement = mean_ss
                #continue

            if key in res_cons_dict.keys():
                cons = res_cons_dict[key]
            else:
                cons = mean_cons

            if key in rsa_agreement_res_dict:
                rsa_agreement = rsa_agreement_res_dict[key]
            else:
                rsa_agreement = mean_rsa
                #continue

            if key in distance_res_dict:
                dist = distance_res_dict[key]
            else:
                dist = mean_dist
                # If the distance value doesn't exist for this residue then skip it
                #continue


            #feature_list.append([target+'~'+model,res_id,dist,cons,ss_agreement,rsa_agreement,\
                             #overall_rsa,overall_ss])



            if train_flag == False:
                test_set.append([target+'~'+model,res_id,dist,cons,ss_agreement,rsa_agreement,                             overall_rsa,overall_ss])
            else:

                train_set.append([target+'~'+model,res_id,dist,cons,ss_agreement,
                                  rsa_agreement,\
                             overall_rsa,overall_ss])
                
    
    
    #print target+'~'+model
    #print 'feature_list=',len(feature_list),'ss_len=',len(ss_agreement_res_dict),\
    #'rsa_len=',len(rsa_agreement_res_dict),'cons_len=',len(res_cons_dict),'dist_len=',len(distance_res_dict)




    # order by the 1st then 2nd column
    test_set = sorted(test_set,key=lambda x: (x[0],x[1])) 
    train_set = sorted(train_set,key=lambda x: (x[0],x[1]))


    # write to fil
    os.chdir('/net/kihara/ialqasse/Desktop/output/data_file')

    #out_file = open('feature.table', 'w')
    test_file = open('test_May4','w')
    train_file = open('train_May4','w')


    #out_file.write('%s %i %f %f %f %f %f %f\n'\
                  #%(target+'~'+model,res_id,dist,cons,ss_agreement,rsa_agreement,\
                  # overall_rsa,overall_ss))


    for item in train_set:
        target_model = item[0]
        res_id = item[1]
        dist = item[2]
        cons = item[3]
        ss_agreement = item[4] 
        rsa_agreement = item[5]
        overall_rsa = item[6]
        overall_ss = item[7]
        train_file.write('%s %i %f %f %f %f %f %f\n'                   %(target_model,res_id,dist,cons,ss_agreement,rsa_agreement,                     overall_rsa,overall_ss))

    for item in test_set:
        target_model = item[0]
        res_id = item[1]
        dist = item[2]
        cons = item[3]
        ss_agreement = item[4] 
        rsa_agreement = item[5]
        overall_rsa = item[6]
        overall_ss = item[7]
        test_file.write('%s %i %f %f %f %f %f %f\n'                       %(target_model,res_id,dist,cons,ss_agreement,rsa_agreement,                         overall_rsa,overall_ss))

    #out_file.close()
    train_file.close()
    test_file.close() 



   


# In[20]:

print ">>>>>>>>>>>>>>>>>count",count

