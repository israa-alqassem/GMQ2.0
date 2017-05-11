
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""auxiliary_function_and_variables.py: 
        contains supplemental functions and variables"""

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, KiharaLab"


# In[15]:

import os
from Bio import pairwise2
fasta_files_dir = '/net/kihara/ialqasse/Desktop/input_files/target_seq'
pdb_files_dir = '/net/kihara/ialqasse/Desktop/input_files/target_pdb'


# In[16]:

# key: 3 letter code, val: single letter
aminoacid_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


# In[17]:

# See STRUDE readme doc
# Mapping between secondary structure (SS) classes in STRIDE and SS in PSIPRED
#H	      AlphaHelix
#G	      310Helix
#I	      PI-helix
#E	      Extended conformation
#B or b   Isolated bridge
#T	      Turn
#C	      Coil (none of the above)

secondary_classes_dict = {'H':'H', 'G': 'H', 'I' : 'H', 
                              'E':'E', 'B':'E', 'b':'E','T' : 'E', 
                              'C': 'C'}


def get_pdb_model_list():
    DIR = '/net/kihara/ialqasse/Desktop/input_files/target_pdb'
    extension = '.mod' # same as '.pdb'
    # the immediate directories
    folder_list = [folder for folder in next(os.walk(DIR))[1]]
    for folder in folder_list:
        if folder == 'modified_pdb':
            continue

        modellist = os.listdir(DIR+'/'+folder)
                
        modellist = [modellist[i].split('.')[0] for i in range(0,len(modellist))]
                
        return modellist
    
    
    

def get_fasta_file_paths():
    """
        This function returns a dictionary of targets as keys and the path to theit fasta files as values
    """
    filelist = []
    extension = '.seq'
    target_fastafile_dict = {}
    for f in os.listdir(fasta_files_dir):
        if f.endswith(extension):
            target_id = f.split('.')[0]
            target_fastafile_dict[target_id] = fasta_files_dir+"/"+f
 
    return target_fastafile_dict


# In[19]:

def get_pdb_file_paths():
    """
        This function returns a dictionary of targets as keys and the path to theit pdb files as values
    """
    DIR = pdb_files_dir
    extension = '.mod' # same as '.pdb'
    # the immediate directories
    folder_list = [folder for folder in next(os.walk(DIR))[1]]
    target_pdbfile_dict = {}
    for folder in folder_list:
        if folder == 'modified_pdb':
            continue
        for f in os.listdir(DIR+'/'+folder):
            if f.endswith(extension) and f.startswith('3D-JIGSAW_V4-0_'):
                target_pdbfile_dict[folder] = DIR+'/'+folder+'/'+f

    return target_pdbfile_dict


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



# In[20]:

def get_fasta_pdb_path_dict():
    """
        This function returns a dictionary of targets as keys and the path to theit fasta and pdb files as values
    """
    target_fastafile_dict = get_fasta_file_paths()
    target_pdbfile_dict = get_pdb_file_paths()
    target_file_dict = {}
    for target in target_pdbfile_dict.keys():
        target_file_dict[target] = [target_fastafile_dict[target],target_pdbfile_dict[target]]
    return target_file_dict


# In[21]:

def diff_letters(str1,str2):
    """
        for two aligned strings str1, str2, this fucntion returns the mismatch count
    """
    str_len = len(str1) if len(str1) <= len(str2) else len(str2)
    return sum ( str1[i] != str2[i] for i in range(str_len))


# In[22]:

def seq_match_accuracy(str1, seq1,str2,seq2):
    """
        parameters:
            str1 and str2: are either surface area strings or secondary classes string (the actual and predicted one respectively)
            seq1 and seq2: are the actual and predicted amino acid sequences
        returns:
        new_str1, new_str2: the same as input string but aligned (dashes used to replace missing or mismatched values)
        seq1: the actual amindo acid sequence 
            
    """
    alignments = pairwise2.align.globalxx(seq1,seq2)
    seq1 =  alignments[0][0]
    seq2 =  alignments[0][1]
    new_str1 = align_output_str_with_residue_seq(str1,seq1)
    new_str2 = align_output_str_with_residue_seq(str2,seq2)
    diff_sum = diff_letters(new_str1,new_str2)
    seq_len = len(seq1)
    return (seq_len - diff_sum)/float(seq_len), new_str1, new_str2,seq1


# In[23]:

def align_output_str_with_residue_seq(str1,seq1):
    """ 
       This function is called from  seq_match_accuracy to align the output string (SA or SS) with amino acid sequence
    """
    new_str1 = ''
    index = 0
    for char in seq1:
        if char == '-':
            new_str1 += '-'
        else:
            new_str1 += str1[index]
            index += 1
    return new_str1



def file_exist(filename):
    if not os.path.exists(filename):
        return False

    if os.stat(filename).st_size == 0:
        return False
    
    return True

