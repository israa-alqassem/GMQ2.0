
# coding: utf-8

# In[1]:

import os
from pandas import *
from Bio import pairwise2
from auxiliary_function_and_variables import *
from RSA_agreement2 import *


# In[2]:

output_per_model_dir = '/net/kihara/ialqasse/Desktop/output/pdb_models_output'
RSAresults_dir = '/net/kihara/ialqasse/Desktop/output/RSA_results/'
model_list = get_pdb_model_list()


# In[3]:

# Run one time to create directories for each model
os.chdir(RSAresults_dir)
for model in model_list:
    os.system('mkdir %s' % model)


# In[4]:

# This function calculates RSA agreement per residue over a window of size 21
def calculate_accuracy_over_window2(out_DIR,target,actualSA_str,predictedSA_str,residue_list,orig_residue_ids):


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
        
    os.chdir(out_DIR)
    filename = target+'_RSA_per_residue.txt'
    header = ['NUM','RES', 'RSA_Agreement']
    f = open(filename, 'w')
    f.write("%s\t%s\t%s\n" % (header[0],header[1],header[2]))
    for item in output_accuracy_list:
        f.write("%d\t%c\t%f\n" % (item[0],item[1],item[2]))
    f.close()


# In[5]:

#if __name__ == '__main__':

target_accfile_dict = get_acc_outputfile_files(accpro_out_output)


for model in model_list:   
    
    #target_rsafile_dict = get_naccess_output_files(naccess_out_dir)
    #already got output for these
    if model.startswith('3D-JIGSAW_'):
        continue
    
 
    print 'Currently working on model ---> ', model
    
    SA_accuracy_output_dict = {}
    
    for target in target_accfile_dict.keys():
        
        #print '\ttarget-->'+target

        #naccess_rsa_file = target_rsafile_dict[target]
        naccess_rsa_file = output_per_model_dir+'/'+model+'/'+target+'/'+model+'.rsa'  
        ACCPro_acc_file = target_accfile_dict[target]
        
        
        if not os.path.exists(naccess_rsa_file):
            print '\tWARNING: pdb file does not exist for: '
            print '\t'+naccess_rsa_file
            continue
            
        if os.stat(naccess_rsa_file).st_size == 0:
            print '\tWARNING: encountered empty file! '
            print '\t'+naccess_rsa_file
            continue
        
        

        actual_SA_str, actual_residue_seq,orig_residue_ids = get_actual_surface_area('',naccess_rsa_file)
        predicted_SA_str,pred_residue_seq = get_predicted_surface_area(target,'',ACCPro_acc_file)

        SA_accuracy_output_dict[target],actual_SA_str2,predicted_SA_str2,actual_residue_seq_modified = seq_match_accuracy(actual_SA_str, actual_residue_seq,predicted_SA_str,pred_residue_seq)


        calculate_accuracy_over_window2(RSAresults_dir+model,target,actual_SA_str2,predicted_SA_str2,actual_residue_seq_modified,orig_residue_ids)



        #i += 1
        #if i==4:
            #break

   
    os.chdir(RSAresults_dir+model)
    filename = 'overall_RSA.txt'
    header = ['Target', 'RSA_Agreement']
    f = open(filename, 'w')
    f.write("%s\t%s\n" % (header[0],header[1]))
    for target in sorted(SA_accuracy_output_dict.keys()):
        f.write("%s\t%s\n" % (target,SA_accuracy_output_dict[target]))

    f.close()
    
    #break

