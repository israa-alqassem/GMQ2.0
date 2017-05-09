
# coding: utf-8

# In[1]:

#!/usr/bin/env python

"""prepare_gmq_input_file.py: 
            prepares GMQ output file (either training or test file depending on the input file)
             the code here is an exact copy of prepare_gmq_training_file.py, I just seperated them to generate training and testing data simultaneously 
            reads feature_table which lists;
            target~model indetifier residueID distance feature1 feature2 etc.
            and prepares GMQ input, run this code to generate both training and testing sets
            
            In this modified version of the code:
            1. I delt with residue IDs all the time (as they appear in the pdb file) 
                the previous code considers neighbor IDs and residue indices (i.e, postions in the array)!!
                Which could be a probolematic if the IDs and indices don't match as in the case of missing residues
            2. I fixed the bug in the previous code, now I only consider the nearest neighbors, previous code did not
            3. here CC case (when the residue and its neighbor both have class c) is given a weight, but in the previous code this case is ignored
            """

__authos__      = "Xuejiao Kang and re-written by Israa Alqassem"
__copyright__   = "Copyright 2017, GMQ"


# In[2]:

import os,sys
out_dir = '/net/kihara/ialqasse/Desktop/output/'


# In[3]:

def get_res_indices_dict(residue_indices):
    '''
     This function maps between residue ids (as they appear in pdb files) to a sequence 
     GMQ assumes that the resiude ids are in consecutive order which may not always be the
     case (when a model has missing residues)
    '''
    res_indices_dict = {}
    for i in range(0, len(residue_indices)):
        res_indices_dict[residue_indices[i]] = i
    return res_indices_dict


# In[4]:

# old func name: name_num_list
def get_targetmodel_residue_count_lists(feature_file):
    '''
        Feature file list the features for each model~target
        This function returns two lists; the number of residues in each model~target
        and all the names of all model~target
    '''
    f = open(feature_file,"r")
    name_residue_count_dict = {} # value: name of model-target, value: residue count
    for row in f:
        row = row.split()
        #print 'row>>>>',row
        if(row[0] not in name_residue_count_dict):
            name_residue_count_dict[row[0]] = 1
        else:
            name_residue_count_dict[row[0]] += 1
            
    f.close()
    
    ## I copied them in lists instead of returning dict.keys() dict.values() because dict doesny
    ## preserve the order of items
    targetmodel_list = []
    residue_count_list = []
    for key in name_residue_count_dict.keys():
        targetmodel_list.append(key)
        residue_count_list.append(name_residue_count_dict[key]) 
        
    return residue_count_list,targetmodel_list


# In[5]:

## old func name: neighbor_feature
def get_nearest_neighbor_residues(neighbor_file,neighbor_cutoff,residue_id,max_neighbor,                         residue_indices_list):
    
    '''
        neighbor_file: it has the extension .ca.dist, it lists the residues with their neighbors and distances
        neighbor_cutoff: neighbor cutoff distance, if a residue within this distance we consider it a neighbor
        residue_id: the id of the residue to find its neighbors
        max_neighbor: max number of neigbour to consider
        residue_indices_list: a list of residue indices read from input feature table
        
        Returns: 
            the residue ids of the nearest neighbouring residues of residue_id
            the ids of the residues that we read from pdb files
    '''

    f = open(neighbor_file,"r")
    neighbor_list = []
    
    for row in f:    
        row = row.split()
        
        if len(row)==9 and row[2]=="CA": # check the format of each line  
            this_res_id = int(row[1])
            neighbor_res_id = int(row[5])
            neighbor_dist = float(row[8])
            
            ########>>>>>>>>>>>>>>>>>>> DOUBLE CHECK THIS
            if this_res_id == residue_id: #residue_indices_list[res_index]:
                if neighbor_dist <= neighbor_cutoff:
                    neighbor_list.append([neighbor_res_id,neighbor_dist])
    f.close()

    #print '**get_nearest_neighbor_residues**'
    #print 'neighbor_list> ',neighbor_list

    # Choose the nearest neighbors if a residue has too many neighbors
    if len(neighbor_list) > max_neighbor:
        neighbor_list.sort(key=lambda x: x[1]) # sort neighbor residues based on distances
        return [neighbor_id[0] for neighbor_id in neighbor_list[0:max_neighbor]]
    else:
        return [neighbor_id[0] for neighbor_id in neighbor_list]
    
    f.close()


# In[6]:

## write_neighbor2file_new
def write_neighbor_to_output_file(sec_struct_filename, neighbor_list, output_f,                            res_id,residue_indices,weight_list):
    
    '''                    
    sec_struct_filename:  a list of res_id, secondary class for a target
    total_num: total number of neighbors
    neighbor_list: a list of neighboring residues ids
    output_f: output file to write to
    res_id: residue id as appeared in pdb file
    residue_indices: residue indices list
    weight_list: has the weights for hh, bb, hb or bh, hc or bc respectively
    
    '''
    
    ## Read secondary structure summary file ##
    sec_struct_file = open(sec_struct_filename,"r")
    
    residue_id_list = []
    residue_class_list = []
    for row in sec_struct_file:
        row = row.split()
        residue_id_list.append(int(row[0]))
        residue_class_list.append(row[1])  # Residue secondary class
        
    sec_struct_file.close()
    
    res_dict = get_res_indices_dict(residue_indices)
    
    for i in range(0,len(neighbor_list)):
        # Check if the neighbor exista in our feature input file
        if neighbor_list[i] in residue_indices: 

            #res_class = residue_class_list[int(residue_id_list.index(residue_indices[res_index]))]
            #neighbor_class = residue_class_list[residue_id_list.index(neighbor_list[i])]
            
            
            ## Get the secondary class of the residue and it's neighbor if both exist in our list
            
            res_index = residue_id_list.index(res_id)
            neighbor_index = residue_id_list.index(neighbor_list[i])
            
            if res_index >= 0 and res_index < len(residue_id_list) and             neighbor_index >= 0 and neighbor_index < len(residue_id_list): 
                res_class = residue_class_list[res_index]
                neighbor_class = residue_class_list[neighbor_index]
            else:
                continue
            
            flag = 0

            if res_class =="H" and neighbor_class =="H":

                weight= weight_list[0] # hh"
                flag=1

            elif res_class=="B" and neighbor_class =="B":

                weight = weight_list[1] # "bb"
                flag=1

            elif (res_class =="H" and neighbor_class =="B") or                       (res_class =="B" and neighbor_class=="H"):

                weight = weight_list[2] # hb or bh
                flag=1

            ## BIG NOTE: the following case considers CC, the orignal code ignores that
            elif res_class == "C" or neighbor_class =="C":

                weight = weight_list[3] # hc, ch, bc, cb, cc
                flag=1

            else:
                print 'Error: Undefined residue or neighbor class: ',res_class,neighbor_class
                
            

            ## BIG NOTE: I considered the residue ids (original ids from pdb file)
            ## old code delt with residue positions/indices in the list
            if flag == 1:
                output_f.write(str(10)+"\n")
                output_f.write(str(2)+" "+str(res_dict[res_id])+" "+str(res_dict[neighbor_list[i]])+"\n")
                output_f.write(str(1)+"\n")
                output_f.write(str(1)+" "+str(1)+" "+str(weight)+"\n")

                output_f.write(str(11)+"\n")      
                output_f.write(str(2)+" "+str(res_dict[res_id])+" "+str(res_dict[neighbor_list[i]])+"\n")
                output_f.write(str(1)+"\n")
                output_f.write(str(0)+" "+str(1)+" "+str(weight)+"\n")

                output_f.write(str(12)+"\n")      
                output_f.write(str(2)+" "+str(res_dict[res_id])+" "+str(res_dict[neighbor_list[i]])+"\n")
                output_f.write(str(1)+"\n")
                output_f.write(str(0)+" "+str(0)+" "+str(weight)+"\n")

                output_f.write(str(13)+"\n")
                output_f.write(str(2)+" "+str(res_dict[res_id])+" "+str(res_dict[neighbor_list[i]])+"\n")
                output_f.write(str(1)+"\n")
                output_f.write(str(1)+" "+str(0)+" "+str(weight)+"\n")



# In[7]:

## old func name: label_file_new
def write_residue_labels_to_output_file(output_f,target_model,total_residue,input_feature,dist_cutoff):
    '''
        This functions writes to the output file a label for each residue
        to identify whether a residue within the specified dist_cutoff or not
            
            output_f: the name of the output file to write to
            target_model: the target~model identifier to write its output to file
            total_residue: total number of residue for this target~model
            input_feature: input feature file to read from
            dist_cutoff: the predefined distance cutoff val
            
        This function also returns a list of residue ids (the ones we read from pdb file)
        not necessarily in consecutive order
            
            
    '''
    input_f = open(input_feature,"r")
    
    #residue_indices = [str(x) for x in range(1,total_residue+1)]
    residue_indices = []
    count = 0
    
    output_f.write(str(total_residue)+'\n')
    
    for row in input_f:
        row = row.split()
        this_target_model = row[0]
        residue_id = int(row[1])   # residue ids as read from pdb file
        distance_val = float(row[2])
        
        ## read data for the specified target model
        if this_target_model == target_model:
            count += 1
            
            residue_indices.append(residue_id)
            
            ## Write labels, i.e., whether this residue within the cutoff distance or not
            if distance_val > dist_cutoff:
                output_f.write('2 0\n')
            else:            
                output_f.write('2 1\n')
      
        if count == total_residue:
            break
            
    input_f.close()
    return residue_indices


# In[8]:

## old func name: feature_file_new
def write_features_to_output_file(output_f,target_model,total_residue,input_feature,neighbor_cutoff_ca,                     num_neighbor,residue_indices,weight_list):
    '''
        output_f: output file
        target_model: target_model name
        total_residue: total number of amino acids for the specified target~model
        input_f: input file
        neighbor_cutoff_ca: neighbor cutoff c_alpha distance
        num_neighbor: max number of neighbors to consider
        residue_indices: list of residue indices for this target~model

    '''
    input_f=open(input_feature,"r")
    
    res_dict = get_res_indices_dict(residue_indices)
    
    res_count = 0
    for row in input_f:
        row = row.split()
        
        if row[0] == target_model:
            # The input feature_table has the format: , 
            # target~model, residue_id, distance, 1st feature, 2nd feature etc.
            res_id = int(row[1])
            fearure_itr = 3 # the feature list starts at the 3rd item in each row
            feature_id = 0
            for j in range(fearure_itr,len(row)):
                output_f.write(str(feature_id)+"\n")
                output_f.write("1 "+str(res_dict[res_id])+"\n")
                output_f.write("1 \n")
                output_f.write("1 "+str(row[j])+"\n")
                feature_id+=1
                output_f.write(str(feature_id)+"\n")
                output_f.write("1 "+str(res_dict[res_id])+"\n")
                output_f.write("1 \n")
                output_f.write("0 "+str(row[j])+"\n")
                feature_id+=1
            

            
            path_structure_ca = out_dir+'/ca_dist'+"/"+row[0][6:] + "/" + row[0][0:5] +".ca.dist"
            path_2ndstructure = out_dir+"/residue_secondary_class_pairs/"+row[0][6:] + "/" + row[0][0:5]+".res.stride"
            
            if os.path.exists(path_structure_ca):
                
                neighbor_list = get_nearest_neighbor_residues(path_structure_ca,neighbor_cutoff_ca,                                                        res_id,num_neighbor,residue_indices)
                
                
                write_neighbor_to_output_file(path_2ndstructure, neighbor_list, output_f,                                            res_id,residue_indices,weight_list)

                
            else:
                print 'Error: '+path_structure_ca+"  does not exist!"
                
                
            res_count+=1
            if res_count==int(total_residue):
                output_f.write("-1\n")
                break
                
    input_f.close()



# In[9]:

if __name__ == '__main__':
    
    #out_dir = '/Users/isra/Desktop/Spring17/KiharaLab/python_code/'
   
    input_feature_file = out_dir+"data_file/test_May4"
    os.chdir(out_dir)
    
    ## TODO:- ask the user to input these features
    ca_cutoff = 5.0        #float(sys.argv[1])
    neighbor_cutoff = 6.0   #float(sys.argv[2])
    num_neighbor = 3        #int(sys.argv[3])
    weight_list = [0.8, 0.9, 1.0, 0.1] ## hh, bb, hb or bh, *c or c* respectively
    
    
    ## get a list of all the unique targets and the totdal num of residues each target has 
    [total_residue_count,target_model_list]=get_targetmodel_residue_count_lists(input_feature_file)

    #print "ca_cutoff: "+str(ca_cutoff)
    #print "neighbor cutoff:  "+str(neighbor_cutoff)
    #print "weight_list: ",weight_list
    #print "Total number of targets in the feature table = ",len(target_model_list)
    
    
    #output_file = open("test4"+str(ca_cutoff)+"_"+str(neighbor_cutoff)+"_"+str(weight_list[0])+"_"+str(weight_list[1])+"_"+str(weight_list[2])+"_"+str(weight_list[3]),"w")    
    output_file = open('gmq_test_may4','w')
    
    count = 0

    for i in range(len(target_model_list)):
        ## Proteins with long sequences of amino acid cause problems in the inference step
        ## reduce their neighbors to 1
        if total_residue_count[i] > 300:
            continue
        else:
            count += 1
            
        print target_model_list[i]

        
        residue_indices = write_residue_labels_to_output_file(output_file,target_model_list[i],
                                                              total_residue_count[i],input_feature_file,
                                                              ca_cutoff)

        write_features_to_output_file(output_file,target_model_list[i],total_residue_count[i],
                                      input_feature_file,neighbor_cutoff,num_neighbor,residue_indices,
                                      weight_list)


    output_file.close()



# In[10]:

print ">>>>>>>",count

