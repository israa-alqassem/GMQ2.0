
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""RSA_agreement.py: 
        calculates relative surface area agreement over a window of 21 residues and over the entire model"""

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, KiharaLab"


# In[1]:

import os
from pandas import *
from get_fasta_pdb_dict import *
from Bio import pairwise2
from auxiliary_function_and_variables import *


# In[2]:

naccess_out_dir = '/net/kihara/ialqasse/Desktop/output/NACCESS_output'
accpro_out_output = '/net/kihara/ialqasse/Desktop/output/ACCPro_output'
out_dir_RSA_per_residue = '/net/kihara/ialqasse/Desktop/output/RSA_per_residue'
out_dir_RSA_all = '/net/kihara/ialqasse/Desktop/output/RSA_all'
target_fastafile_dict = get_fasta_file_paths()


# In[3]:

def get_residue_seq_from_fastafile(filename):
    residue_seq = open(filename, 'r').read().split('\n')[1]
    return residue_seq


# In[4]:

def get_naccess_output_files(DIR,extension = '.rsa'):
    # the immediate directories
    folder_list = [folder for folder in next(os.walk(DIR))[1]]
    target_naccess_rsa_dict = {}
    for folder in folder_list:
        for f in os.listdir(DIR+'/'+folder):
            if f.endswith(extension): #and f.startswith('3D-JIGSAW_V4-0_'):
                target_naccess_rsa_dict[folder] = DIR+'/'+folder+'/'+f

    return target_naccess_rsa_dict


# In[5]:

def get_acc_outputfile_files(DIR,extension = '.acc'):
    target_accfile_dict = {}
    for f in os.listdir(DIR):
        if f.endswith(extension):
            target_id = f.split('.')[0]
            target_accfile_dict[target_id] = DIR+"/"+f

    return target_accfile_dict


# In[6]:

def get_actual_surface_area(DIR,filename):
    naccess_output = read_csv(DIR+filename, delim_whitespace=True, skiprows=[0,1,2,3],header=None)
    actual_surface_area_str = ''
    residue_id_list = []
    residue_seq = ''
    count=0
    for _,row in naccess_output.iterrows():
        if row[0] == 'RES':
            count+=1
            residue_name = aminoacid_dict[row[1]]
            residue_num = int(row[2])
            residue_id_list.append(residue_num)
            residue_seq += residue_name
            relative_accessibility = float(row[4])
            
            if relative_accessibility >= 25.:
                actual_surface_area_str += 'e' #expected
            else:
                actual_surface_area_str += '-'
           
    return actual_surface_area_str,residue_seq,residue_id_list


# In[7]:

def get_predicted_surface_area(target,DIR,filename):
    accpro_output = open(DIR+ACCPro_acc_file, 'r').read()
    #print len(accpro_output.split('\n')[1])
    predicted_surface_area_str = accpro_output.split('\n')[1]
    
    fasta_file = target_fastafile_dict[target]
    residue_seq = get_residue_seq_from_fastafile(fasta_file)
        
    return predicted_surface_area_str,residue_seq


# In[8]:

# This function calculates RSA agreement per residue over a window of size 21
def calculate_accuracy_over_window(target,actualSA_str,predictedSA_str,residue_list,orig_residue_ids):


    output_accuracy_list = []
    total_residues = len(residue_list)
    i=0         # residue index
    count = 0
    for residue in residue_list:
        #print 'i = ',i
        if i <=10 and i+10 <= total_residues: #residue located at the begining, window_size < 21
            start = 0
            end = i+10

        elif i > 10 and i+10 <= total_residues: #residue located in the middle, window_size = 21
            start = i-10
            end = i+10

        elif i+10 < total_residues: #residue located at the end, window_size < 21
            start = i-10
            end = total_residues

        diff = diff_letters(actualSA_str[start:end],predictedSA_str[start:end])
        win_size = len(actualSA_str[start:end])+1

        if residue != '-':
            output_accuracy_list.append([orig_residue_ids[count],residue,(win_size-diff)/float(win_size)])
            count += 1

        i += 1


    #return output_accuracy_list
        
    os.chdir(out_dir_RSA_per_residue)
    filename = target+'_RSA_per_residue.txt'
    header = ['NUM','RES', 'RSA_Agreement_over_window_of_21_residues']
    f = open(filename, 'w')
    f.write("%s\t%s\t%s\n" % (header[0],header[1],header[2]))
    for item in output_accuracy_list:
        f.write("%d\t%c\t%f\n" % (item[0],item[1],item[2]))
    f.close()


# In[9]:

if __name__ == '__main__':
   
    target_rsafile_dict = get_naccess_output_files(naccess_out_dir)
    target_accfile_dict = get_acc_outputfile_files(accpro_out_output)
    SA_accuracy_output_dict = {}
    i = 0
    for target in target_rsafile_dict.keys():
        naccess_rsa_file = target_rsafile_dict[target]
        ACCPro_acc_file = target_accfile_dict[target]

        actual_SA_str, actual_residue_seq,orig_residue_ids = get_actual_surface_area('',naccess_rsa_file)
        predicted_SA_str,pred_residue_seq = get_predicted_surface_area(target,'',ACCPro_acc_file)

        SA_accuracy_output_dict[target],actual_SA_str2,predicted_SA_str2,actual_residue_seq_modified = seq_match_accuracy(actual_SA_str, actual_residue_seq,predicted_SA_str,pred_residue_seq)


        calculate_accuracy_over_window(target,actual_SA_str2,predicted_SA_str2,actual_residue_seq_modified,orig_residue_ids)



        #i += 1
        #if i==4:
            #break


    os.chdir(out_dir_RSA_all)
    filename = 'overall_RSA.txt'
    header = ['Target', 'SS_Agreement']
    f = open(filename, 'w')
    f.write("%s\t%s\n" % (header[0],header[1]))
    for target in sorted(SA_accuracy_output_dict.keys()):
        f.write("%s\t%s\n" % (target,SA_accuracy_output_dict[target]))

    f.close()


