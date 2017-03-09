
# coding: utf-8

# In[1]:

import os
from pandas import *
#from get_fasta_pdb_dict import *
from Bio import pairwise2
from auxiliary_function_and_variables import *


# In[2]:

stride_out_dir = '/net/kihara/ialqasse/Desktop/output/STRIDE_output'
psipred_out_dir = '/net/kihara/ialqasse/Desktop/output/PSIPRED_output'
out_dir_SS_per_residue = '/net/kihara/ialqasse/Desktop/output/SS_per_residue'
out_dir_SS_all = '/net/kihara/ialqasse/Desktop/output/SS_all'


# In[3]:

def get_stride_output_paths(DIR,extension = ''):
    DIR = stride_out_dir
    target_strideout_dict = {}
    for f in os.listdir(DIR):
        target_strideout_dict[f] = DIR+"/"+f
        

    return target_strideout_dict


# In[4]:

def get_psipred_output_paths(DIR,extension = '.ss'):
    target_psipredout_dict = {}
    for f in os.listdir(DIR):
        if f.endswith(extension):
            target_id = f.split('.')[0]
            target_psipredout_dict[target_id] = DIR+"/"+f


    return target_psipredout_dict


# In[5]:

def read_stride_data(DIR,filename):
    start_reading = False
    stride_ss_output_dict = {}
    residue_seq = ''
    residue_id_list = []
    secondary_class_str = ''
    
    with open(DIR+stride_file) as f:
        for line in f:
            if 'Detailed secondary structure assignment' in line:
                start_reading = True

            if start_reading == True:
                if line[0:3] == 'ASG':
                    line = line.split()  
                    stride_id,stride_aminoacid,stride_secondary_class = line[3],aminoacid_dict[line[1]],secondary_classes_dict[line[5]]
                    stride_ss_output_dict[stride_id+stride_aminoacid] = stride_secondary_class
                    residue_seq += stride_aminoacid
                    residue_id_list.append(int(stride_id))
                    secondary_class_str += stride_secondary_class
                    

    return stride_ss_output_dict,residue_seq,residue_id_list,secondary_class_str


# In[6]:

def read_psipred_data(DIR, filename):
    psipred_data= read_csv(DIR+filename,delim_whitespace=True,header=None)
    psipred_ss_output_dict = {}
    residue_seq = ''
    residue_id_list = []
    secondary_class_str = ''
    for i in range(len(psipred_data)):     
        psipred_id,psipred_aminoacid,psipred_secondary_class = str(psipred_data[0][i]),psipred_data[1][i],psipred_data[2][i]
        psipred_ss_output_dict[psipred_id+psipred_aminoacid] = psipred_secondary_class
        
        residue_seq += psipred_aminoacid
        residue_id_list.append(psipred_id)
        secondary_class_str += psipred_secondary_class
        
    return psipred_ss_output_dict,residue_seq,residue_id_list, secondary_class_str


# In[16]:

# This function calculates secondary structure agreement per residue over a window of size 21
def calculate_accuracy_over_window(target,actualSS_str,predictedSS_str,residue_list,orig_residue_ids):


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
    os.chdir(out_dir_SS_per_residue)
    filename = target+'_SS_per_residue.txt'
    header = ['NUM','RES','Real_class', 'Predicted_class' ,'SS_agreement']
    f = open(filename, 'w')
    f.write("%s\t%s\t%s\t%s\t%s\n" % (header[0],header[1],header[2],header[3],header[4]))
    for item in output_accuracy_list:
        f.write("%d\t%c\t%c\t%c\t%f\n" % (item[0],item[1],item[2],item[3],item[4]))
    f.close()


# In[17]:

if __name__ == '__main__':
    
    target_stride_out_dict = get_stride_output_paths(stride_out_dir)
    target_psipred_out_dict = get_psipred_output_paths(psipred_out_dir)
    SS_output_dict = {}
    
    for target in target_stride_out_dict.keys():

        stride_file = target_stride_out_dict[target]
        psipred_ss_file = target_psipred_out_dict[target]

        stride_ss_output_dict,stride_residue_seq,stride_residue_id_list, stride_SS_str = read_stride_data('',stride_file)
        psipred_ss_output_dict,psipred_residue_seq,psipred_residue_id_list, psipred_SS_str = read_psipred_data('',psipred_ss_file)

        SS_output_dict[target],actual_SS_str2,predicted_SS_str2,actual_residue_seq_modified = seq_match_accuracy(stride_SS_str, stride_residue_seq,psipred_SS_str,psipred_residue_seq)


        output_accuracy_list = calculate_accuracy_over_window(target,actual_SS_str2,predicted_SS_str2,actual_residue_seq_modified,stride_residue_id_list)
        
        


    #os.chdir(out_dir_SS_all)
    #filename = 'overall_SS.txt'
    #header = ['Target', 'SS_Agreement']
    #f = open(filename, 'w')
    #f.write("%s\t%s\n" % (header[0],header[1]))
    #for target in sorted(SS_output_dict.keys()):
        #f.write("%s\t%s\n" % (target,SS_output_dict[target]))

    #f.close()


# In[ ]:



