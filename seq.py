# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 19:05:42 2021

@author: Chaewon
"""

#%% read TF-seq csv files and append to a pandas DataFrame
"""
TFseq_df =
                              Name      Seq
        0      CACCC-bindingfactor   GGGTGG
        1                      ETF   CCACCC
        2                gammaCAC1    CACCC
        3                gammaCAC2    CACCC
        ...
# of rows : 19911
"""

import pandas as pd

TFseq_df = pd.DataFrame()
for i in range(1,20):
    file_location = r'C:\Users\Chaewon\Desktop\TF data\For database'
    file_name = str(i)+r'_sorted.csv'
    path = file_location + r'\\' + file_name
    raw_data = pd.read_csv(path)
    
    temp_df = pd.DataFrame(raw_data,columns=['Name',' Seq'])
    TFseq_df = TFseq_df.append(temp_df, ignore_index = True)

print(TFseq_df)

file_location = r'C:\Users\Chaewon\Desktop\TF data\For database'
file_name = r'SIRT3-sorted.csv'
path = file_location + r'\\' + file_name
raw_data = pd.read_csv(path)

temp_df = pd.DataFrame(raw_data,columns=['Name',' Seq'])
TFseq_df = TFseq_df.append(temp_df, ignore_index = True)

#%% done
"""
seq_list = ['AAAA','C','GG',... ], length : 735
"""
seq_set = set()
for index, row in TFseq_df.iterrows():
    seq_set.add(row[' Seq'])
seq_list = sorted(list(seq_set))

#%% done
"""
seq_TF_dict = {'AAAC':['Sp1','Sp2',... ],
               ...
              'CCC':['Sp4']}
"""
seq_TF_dict = {}
for seq in seq_set: # seq = 'AAAC'
    matched_TFseq_df = TFseq_df[(TFseq_df[' Seq'] == seq)]
    # matched_TF_list = ['Sp1','Sp2',... ]
    matched_TF_set = set()
    for index, row in matched_TFseq_df.iterrows():
        matched_TF_set.add(row['Name'])
    matched_TF_list = sorted(list(matched_TF_set))
    # append {'AAAC': ['Sp1','Sp2',... ]}
    seq_TF_dict[seq] = matched_TF_list
    
#%% done
"""
seq_list2 = ['AAACA','AAGAT','AATCT',... ]
"""    
seq_list2 = sorted(seq_list, key=len)

#%% done
"""
len_list = [5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 19, 20, 22, 36, 48]
len_seq_list_dict = {5: ['AAACA','AAGAT','AATCT',... ],6: []}
"""
len_set = set()
for seq in seq_list2:
    len_set.add(len(seq))
len_list = list(len_set)

len_seq_list_dict = {}
# initialize
for l in len_set:
    len_seq_list_dict[l] = []
# 
for seq in seq_list2:
    len_seq_list_dict[len(seq)].append(seq)
    
#%% 
"""
len_list = [5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 19, 20, 22, 36, 48]
len_seq_list_dict = {5: ['AAACA','AAGAT','AATCT',... ],6: []}
seq_list2 = ['AAACA','AAGAT','AATCT',... ]
seq_TF_dict = {'AAAC':['Sp1','Sp2',... ],
               ...
              'CCC':['Sp4']}

state_rule_dict = {'':['',{'A':'A','C':'C','G':'G','T':'T'}],
                   'A':['',{}],
                   'C':['',{}],
                   'AC':['C',{'A':0,'C':0,'G':0,'T':0}],
                    ...
                   'AAAC':['AC',{'A':0,'C':0,'G':3,'T':0}]}        
"""

"""
pseudo-code inside while loop (tentative)
- pop queue
- 'A','C','G','T' insert 'X' in queue and state_list if 'X' is prefix of any seq in seq_list2


queue = ['A','C','G','T']
state_list = ['']
state_info_dict = {'':['',{'A':'A','C':'C','G':'G','T':'T'}],
                   'A':['',{'A':'A','C':}}]}
"""
def is_prefix(state_cand, len_list, len_seq_list_dict):
    for num in len_list:
        state_length = len(state_cand)
        if state_length <= num:
            num_length_seq_list = len_seq_list_dict[num]
            for seq in num_length_seq_list:
                if state_cand == seq[:state_length]:
                    return 1
    return 0

def is_TFseq(seq, seq_TF_dict):
    try:
        seq_TF_dict[seq]
        return 1
    except:
        return 0
    
def find_child(seq, seq_TF_dict):
    for n in range(1,len(seq)):
        if is_TFseq(seq[n:], seq_TF_dict) == 1:
            return seq[n:]
    return ''

def find_state(seq, state_info_dict):
    for n in range(len(seq)):
        try:
            state_info_dict[seq[n:]]
            return seq[n:]
        except:
            pass
    return ''

state_info_dict={}
state_to_add_in_dict_queue = []
state_to_add_in_dict_queue.insert(0, '')
state_list = []
state_list.append('')

while len(state_to_add_in_dict_queue):
    # add next states in state_list and state_to_add_in_dict_queue
    state = state_to_add_in_dict_queue.pop()
    for base in ['A','C','G','T']:
        next_state_cand = state + base
        if is_prefix(next_state_cand, len_list, len_seq_list_dict):
            state_to_add_in_dict_queue.insert(0, next_state_cand)
            state_list.append(next_state_cand)
            state_info_dict[next_state_cand] = []
    # add {state: info} in state_info_dict
    child = find_child(state, seq_TF_dict)
    transition_dict = {}
    for base in ['A','C','G','T']:
        next_state_cand = state + base
        transition_dict[base] = find_state(next_state_cand, state_info_dict)
    state_info_dict[state] = [child, transition_dict]
    
#%%    
file = open(r'C:\Users\Chaewon\Desktop\TF data\For database\SIRT3_promoter.txt', "r")
promoter = file.read()
promoter_rev_list = []

for ch in promoter:
    if ch == 'A':
        promoter_rev_list.insert(0,'T')
    elif ch == 'C':
        promoter_rev_list.insert(0,'G')        
    elif ch == 'G':
        promoter_rev_list.insert(0,'C')        
    elif ch == 'T':
        promoter_rev_list.insert(0,'A')

promoter_rev = ''.join(promoter_rev_list)
        
#%% 
import pandas as pd
match_result_df = pd.DataFrame(columns=['Name',' Pos',' Seq'])

state = ''
pos = 0

for ch in promoter:
    pos = pos + 1
    state = state_info_dict[state][1][ch]
    #print(pos, ch, state)
    try:
        for TF in seq_TF_dict[state]:
            data = [TF, str(pos)+' (+)', state]
            match_result_df.loc[len(match_result_df.index)] = data
    except:
        pass
    
    child = state_info_dict[state][0]
    while child != '':
        for TF in seq_TF_dict[child]:
            data = [TF, str(pos)+' (+)', child]
            match_result_df.loc[len(match_result_df.index)] = data
        child = state_info_dict[child][0]
        
#%%
state = ''
pos = 0

for ch in promoter_rev:
    pos = pos + 1
    state = state_info_dict[state][1][ch]
    try:
        for TF in seq_TF_dict[state]:
            data = [TF, str(len(promoter)-pos+1)+' (-)', state]
            match_result_df.loc[len(match_result_df.index)] = data
    except:
        pass
    
    child = state_info_dict[state][0]
    while child != '':
        for TF in seq_TF_dict[child]:
            data = [TF, str(len(promoter)-pos+1)+' (-)', child]
            match_result_df.loc[len(match_result_df.index)] = data
        child = state_info_dict[child][0]
        
#%%
match_result_df.to_csv(r'C:\Users\Chaewon\Desktop\TF data\For database\result.csv')

