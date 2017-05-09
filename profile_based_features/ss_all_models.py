
# coding: utf-8

# In[1]:

import os
from pandas import *
from Bio import pairwise2
from auxiliary_function_and_variables import *
from SS_agreement2 import *


# In[2]:

output_per_model_dir = '/net/kihara/ialqasse/Desktop/output/pdb_models_output'
SSresults_dir = '/net/kihara/ialqasse/Desktop/output/SS_results/'
model_list = get_pdb_model_list()


# In[3]:

# Run one time to create directories for each model
os.chdir(SSresults_dir)
for model in model_list:           
        os.system('mkdir %s' % model)


# In[4]:

# This function calculates secondary structure agreement per residue over a window of size 21
def calculate_accuracy_over_window2(out_DIR, target,actualSS_str,predictedSS_str,residue_list,orig_residue_ids):


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

        diff = diff_letters(actualSS_str[start:end],predictedSS_str[start:end])
        win_size = len(actualSS_str[start:end])+1
        
        actual_class = actualSS_str[i]
        predicted_class = predictedSS_str[i]

        if residue != '-':
            output_accuracy_list.append([orig_residue_ids[count],residue,actual_class,predicted_class,
                                         (win_size-diff)/float(win_size)])
            count += 1

        i += 1


    #return output_accuracy_list
    os.chdir(out_DIR)
    filename = target+'_SS_per_residue.txt'
    header = ['NUM','RES','Real_class', 'Pred_class' ,'SS_agreement']
    f = open(filename, 'w')
    f.write("%s\t%s\t%s\t%s\t%s\n" % (header[0],header[1],header[2],header[3],header[4]))
    for item in output_accuracy_list:
        f.write("%d\t%c\t%c\t%c\t%f\n" % (item[0],item[1],item[2],item[3],item[4]))
    f.close()


# In[5]:

#if __name__ == '__main__':

target_psipred_out_dict = get_psipred_output_paths(psipred_out_dir)

for model in model_list:   
    
    #target_stride_out_dict = get_stride_output_paths(stride_out_dir)
    
    if model.startswith('3D-JIGSAW_'):
        continue
        
    print 'Currently working on model ---> ', model
    
    SS_output_dict = {}
    
    for target in target_psipred_out_dict.keys():

        #stride_file = target_stride_out_dict[target]
        stride_file = output_per_model_dir+'/'+model+'/'+target+'/'+target    
        psipred_ss_file = target_psipred_out_dict[target]
        
        
        if not os.path.exists(stride_file):
            print '\tWARNING: pdb file does not exist for: '
            print '\t'+stride_file
            continue
        
        

        stride_ss_output_dict,stride_residue_seq,stride_residue_id_list, stride_SS_str = read_stride_data('',stride_file)
        psipred_ss_output_dict,psipred_residue_seq,psipred_residue_id_list, psipred_SS_str = read_psipred_data('',psipred_ss_file)

        SS_output_dict[target],actual_SS_str2,predicted_SS_str2,actual_residue_seq_modified = seq_match_accuracy(stride_SS_str, stride_residue_seq,psipred_SS_str,psipred_residue_seq)


        output_accuracy_list = calculate_accuracy_over_window2(SSresults_dir+model,target,actual_SS_str2,predicted_SS_str2,actual_residue_seq_modified,stride_residue_id_list)
        
        #break

    
    os.chdir(SSresults_dir+model)
    filename = 'overall_SS.txt'
    header = ['Target', 'SS_Agreement']
    f = open(filename, 'w')
    f.write("%s\t%s\n" % (header[0],header[1]))
    for target in sorted(SS_output_dict.keys()):
        f.write("%s\t%s\n" % (target,SS_output_dict[target]))

    f.close()
    


# In[ ]:



