
# coding: utf-8

# In[ ]:

"""plot_errorbars.py: 
        plots error bars for the data generated on May 4, 2017
        
        
        Training and test data can be found here:
            /net/kihara/ialqasse/Desktop/output/data_file
        
        GMQ output can be found here:
            /net/kihara/ialqasse/Desktop/output/GMQ_output/GMQ_May4
      
Description of the error bars:
The x-axis shows the actual label vs. the predicted label by GMQ (ca_cutoff=5A,neighbor_cutoff=6,num_neighbor=3)
we have 4 classes here [0 vs. 0,0 vs. 1,1 vs. 0,1 vs.1] which correspond to TN,FP,FN,TP, respectively.
while the bars represent the standard deviations, and the little squares are the mean values. 
Error bars help illustrating how well the accuracy function (e.g., RSA, SS, cons) describe the labels, 
i.e., (1) are the secondary structure values higher for the true-positive amd false-postive classes and lower for
true negative and false negative?
(2) are the conservation values smaller for TP class? Since conservation is an entropy we'd expect smaller entropy
for TP class?

P.S. cons values are the values of the 2nd to last column in PSIBlast output
P.S.S get_ss_agreement_per_residue,get_rsa_agreement_per_residue, get_cons_per_residue
all return lists instead of dict because the residues are ordered.

"""

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, KiharaLab"


# In[1]:

get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import sys,os
import numpy
import random

sys.path.append('../')
from auxiliary_function_and_variables import *


# In[2]:

ss_out_dir = '/net/kihara/ialqasse/Desktop/output/SS_results/'
rsa_out_dir = '/net/kihara/ialqasse/Desktop/output/RSA_results/'
psiblast_out_dir ='/net/kihara/ialqasse/Desktop/output/PSIBlast_output/'
calpha_dir = '/net/kihara/ialqasse/Desktop/CALPHA_result'
SS_figures_dir = '/net/kihara/ialqasse/Desktop/output/plots/SS_vs_label/'
gmq_output_file = '/net/kihara/ialqasse/Desktop/output/GMQ_output/GMQ_May4/output'


# In[3]:

def get_ss_agreement_per_residue(model,target):
    
    #ss_agreement_res_list = {}
    ss_list = []
    mean = 0.
    
    filename = ss_out_dir+'/'+model+'/'+target+'_SS_per_residue.txt'
    
    if not file_exist(filename):
        return ss_list,mean

    with open(filename) as f:
        next(f) # ignore header line
        for line in f:
            line = line.split()
            res_id = int(line[0])
            #res = line[1]
            ss_agreement = float(line[4])
            #ss_agreement_res_dict[res_id] = ss_agreement
            ss_list.append(ss_agreement)
            mean += ss_agreement
    mean /= float(len(ss_list))  
    return ss_list, mean


# In[4]:

def get_rsa_agreement_per_residue(model,target):
    
    #rsa_agreement_res_dict = {}
    rsa_agreement_list = []
    mean = 0.
    
    filename = rsa_out_dir+'/'+model+'/'+target+'_RSA_per_residue.txt'
    
    if not file_exist(filename):
        #return rsa_agreement_res_dict,mean
        return rsa_agreement_list,mean

    with open(filename) as f:
        next(f) # ignore header line
        for line in f:
            line = line.split()
            res_id = int(line[0])
            #res = line[1]
            rsa_agreement = float(line[2])
            #rsa_agreement_res_dict[res_id] = rsa_agreement
            rsa_agreement_list.append(rsa_agreement)
            mean += rsa_agreement
            
    mean /= float(len(rsa_agreement_list))
    return rsa_agreement_list, mean


# In[5]:

def get_cons_per_residue(target):
    filename = psiblast_out_dir+target+'_PSSM'
    #cons_res_dict = {}
    cons_list = []
    mean = 0.
    
    if not file_exist(filename):
        #return cons_res_dict,mean
        return cons_list,mean
    
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

                #res = line[1]
                cons = float(line[42])
                #print res_id,
                #print res,
                #print cons
                #cons_res_dict[int(res_id)] = cons
                cons_list.append(cons)
                mean += cons
                
    mean /= float(len(cons_list))
    return cons_list,mean


# In[6]:

def plot_errorbar(feature_values,target_model,actual_label_list,pred_label_list,DIR,fig_title):
    list00 = []
    list01 = []
    list10 = []
    list11 = []

    for i in range(len(feature_values)):
        #print key,val
        val = feature_values[i]
        
        if actual_label_list[i] == 0 and pred_label_list[i] == 0:
            list00.append(val)
        elif actual_label_list[i] == 0 and pred_label_list[i] == 1:
            list01.append(val)
        elif actual_label_list[i] == 1 and pred_label_list[i] == 0:
            list10.append(val)
        elif actual_label_list[i] == 1 and pred_label_list[i] == 1:
            list11.append(val)

                

    mean_vals = numpy.array([numpy.mean(list00),numpy.mean(list01),numpy.mean(list10),numpy.mean(list11)],dtype=float)
    std_vals = numpy.array([numpy.std(list00),numpy.std(list01),numpy.std(list10),numpy.std(list11)],dtype=float)
    labels = numpy.array([1,2,3,4])
    _xticks = ['0 vs. 0','0 vs. 1','1 vs. 0','1 vs. 1']
    plt.xticks(labels, _xticks,fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Actual Label vs. predicted Label',fontsize=12)
    plt.ylabel(fig_title,fontsize=12)
    plt.title(target_model,fontsize=16)

    plt.errorbar(labels,mean_vals,std_vals,linestyle='None',marker='s',ecolor='blue',color='blue')


    plt.savefig(DIR+target_model+".png")
    #plt.show()
    plt.clf()


# In[7]:

test_target_model_list = 'test_models.txt'
target_list = []
f = open(test_target_model_list,'r')
for row in f:
    row = row.split()
    target_list.append(row[0])
f.close()


# In[8]:

SS_figures_dir = '/net/kihara/ialqasse/Desktop/output/plots/SS_vs_label/'
fig_title = 'Secondary Structure Agreement'

actual_label_list = []
pred_label_list = []
target_id = 0

target_model = target_list[target_id].split('~')
feature_list, mean = get_ss_agreement_per_residue(target_model[1],target_model[0])

f = open(gmq_output_file,'r')
for row in f:
    row = row.split()
    len_row = len(row)
    
    if len_row == 1:
        total_res = int(row[0])
        #print 'Target',target_list[target_id],'Total amindo acids =', total_res
    elif len_row  == 2:
        actual_label_list.append(int(row[0]))
        pred_label_list.append(int(row[1]))
        
    elif len_row == 7:
        plot_errorbar(feature_list,target_list[target_id],actual_label_list,pred_label_list,SS_figures_dir,fig_title)
        target_id += 1
        actual_label_list = []
        pred_label_list = []
        feature_list = []
        target_model = target_list[target_id].split('~')
        feature_list, mean = get_ss_agreement_per_residue(target_model[1],target_model[0])

        #print '>>>len(feature_list)=',len(feature_list)
    
    elif len_row == 0:
        continue
    else:
        print ">>>Row of undefined length!", row
        break
        
    if target_id==100:
        break
    
        
f.close()


# In[9]:

RSA_figures_dir = '/net/kihara/ialqasse/Desktop/output/plots/RSA_vs_label/'
fig_title = 'RSA Agreement'


actual_label_list = []
pred_label_list = []
target_id = 0

target_model = target_list[target_id].split('~')
feature_list, mean = get_rsa_agreement_per_residue(target_model[1],target_model[0])


f = open(gmq_output_file,'r')
for row in f:
    row = row.split()
    len_row = len(row)
    
    if len_row == 1:
        total_res = int(row[0])
        #print 'Target',target_list[target_id],'Total amindo acids =', total_res
    elif len_row  == 2:
        actual_label_list.append(int(row[0]))
        pred_label_list.append(int(row[1]))
        
    elif len_row == 7:
        plot_errorbar(feature_list,target_list[target_id],actual_label_list,pred_label_list,RSA_figures_dir,fig_title)
        target_id += 1
        actual_label_list = []
        pred_label_list = []
        feature_list = []
        target_model = target_list[target_id].split('~')
        feature_list, mean = get_rsa_agreement_per_residue(target_model[1],target_model[0])
        
        #print '>>>len(feature_list)=',len(feature_list)
    
    elif len_row == 0:
        continue
    else:
        print ">>>Row of undefined length!", row
        break
        
    if target_id==100:
        break
    
        
f.close()


# In[10]:

CONS_figures_dir = '/net/kihara/ialqasse/Desktop/output/plots/CONS_vs_label/'
fig_title = 'Conservation'

actual_label_list = []
pred_label_list = []
target_id = 0

target_model = target_list[target_id].split('~')
feature_list, mean = get_cons_per_residue(target_model[0])

f = open(gmq_output_file,'r')
for row in f:
    row = row.split()
    len_row = len(row)
    
    if len_row == 1:
        total_res = int(row[0])
        #print 'Target',target_list[target_id],'Total amindo acids =', total_res
    elif len_row  == 2:
        actual_label_list.append(int(row[0]))
        pred_label_list.append(int(row[1]))
        
    elif len_row == 7:
        plot_errorbar(feature_list,target_list[target_id],actual_label_list,pred_label_list,CONS_figures_dir,fig_title)
        target_id += 1
        actual_label_list = []
        pred_label_list = []
        feature_list = []
        target_model = target_list[target_id].split('~')
        feature_list, mean = get_cons_per_residue(target_model[0])
        #print '>>>len(feature_list)=',len(feature_list)
    
    elif len_row == 0:
        continue
    else:
        print ">>>Row of undefined length!", row
        break
        
    if target_id==100:
        break
    
        
f.close()


# In[ ]:



