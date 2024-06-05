#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 21:07:00 2024

@author: mfperez
"""

#%%
import pandas as pd
import decoupler as dc
import os

from pathlib import Path

#%%

home_directory = os.path.expanduser("~")

working_directory = os.path.join(home_directory, "Cel_GRN_manuscript")

os.chdir(working_directory)

output_directory = Path("output/benchmark_out")

output_directory.mkdir(parents = True, exist_ok = True)

#%%

benchRAW = pd.read_csv("output/benchmark_DEstats_RAW.txt",
                  sep = '\t',
                  index_col = 0)

benchRAPToR = pd.read_csv("output/benchmark_DEstats_RAPToR.txt",
                  sep = '\t',
                  index_col = 0)

#%%
obs1 = pd.read_csv("output/benchmark_observations.txt",
                  sep = '\t',
                     index_col=0,
                     encoding='unicode_escape')

#%%

cutoffs_vec = [100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000]


for cutoff in cutoffs_vec:
    
    file_name = "output/GRNs/FIMO_nohomo_" + str(cutoff) + ".txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
    
    thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
    
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            },
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            },
        }
        
    thisoutput = dc.benchmark(benchRAW, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/FIMO_nohomo_cutoff" + str(cutoff) + "_benchRAW.tsv", sep='\t', index=False) 

#%%    

for cutoff in cutoffs_vec:
    
    file_name = "output/GRNs/FIMO_nohomo_" + str(cutoff) + ".txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
    
    thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
    
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            },
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            },
        }
        
    thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
    
    # dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/FIMO_nohomo_cutoff" + str(cutoff) + "_benchRAPToR.tsv", sep='\t', index=False) 

#%%   

cutoffs_vec = [100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000]

for cutoff in cutoffs_vec:

    file_name = "output/GRNs/FIMO_nohomo_" + str(cutoff) + ".txt"
    
    thisGRN_dict = {}
    
    thisGRN_kws = {
        
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            },
        }
        
    shuffle_stats = []
    
    for i in range(1, 101):

        thisGRN_dict[str(cutoff)] = pd.read_table(file_name) 
    
        thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
        
        shuffledict = {str(cutoff) + "_shuffle": thisGRN_dict[str(cutoff) + "_shuffle"]}
        
        thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws) 

        my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
        shuffle_stats.append(my_selection)

    combined_df = pd.concat(shuffle_stats, ignore_index = True)

    combined_df.to_csv("output/benchmark_out/FIMO_nohomo_cutoff" + str(cutoff) + "_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 
    


#%%   

homoparam_vec = [1, 3, 5, 7]

for homoparam in homoparam_vec:
    
    for cutoff in cutoffs_vec:
        
        file_name = "output/GRNs/FIMO_homoparam" + str(homoparam) + "_" + str(cutoff) + ".txt"
        
        thisGRN_dict = {}
        
        thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
        
        thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
        
        thisGRN_kws = {
            
            str(cutoff):{
                'methods' : ['mlm', 'ulm', 'wsum'],
                'consensus' : True,
                # 'dense' : True
                },
            str(cutoff) + "_shuffle":{
                'methods' : ['mlm', 'ulm', 'wsum'],
                'consensus' : True,
                # 'dense' : True
                },
            }
        
        thisoutput = dc.benchmark(benchRAW, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws) 
        
        thisoutput.to_csv("output/benchmark_out/FIMO_homoparam" +str(homoparam) + "_cutoff" + str(cutoff) + "_benchRAW.tsv", sep='\t', index=False) 



#%%   

homoparam_vec = [1, 3, 5, 7]

for homoparam in homoparam_vec:
    
    for cutoff in cutoffs_vec:
        
        file_name = "output/GRNs/FIMO_homoparam" + str(homoparam) + "_" + str(cutoff) + ".txt"
        
        thisGRN_dict = {}
        
        thisGRN_dict[str(cutoff)] = pd.read_table(file_name)
        
        thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
        
        thisGRN_kws = {
            
            str(cutoff):{
                'methods' : ['mlm', 'ulm', 'wsum'],
                'consensus' : True,
                # 'dense' : True
                },
            str(cutoff) + "_shuffle":{
                'methods' : ['mlm', 'ulm', 'wsum'],
                'consensus' : True,
                # 'dense' : True
                },
            }
        
        thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws) 
        
        thisoutput.to_csv("output/benchmark_out/FIMO_homoparam" +str(homoparam) + "_cutoff" + str(cutoff) + "_benchRAPToR.tsv", sep='\t', index=False) 


#%%   

homo_cutoffs_vec = [100, 500, 1000, 2000]

for cutoff in homo_cutoffs_vec:

    file_name = "output/GRNs/FIMO_homoparam5_" + str(cutoff) + ".txt"
    
    thisGRN_dict = {}
    
    thisGRN_kws = {
        
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            },
        }
        
    shuffle_stats = []
    
    for i in range(1, 101):

        thisGRN_dict[str(cutoff)] = pd.read_table(file_name) 
    
        thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
        
        shuffledict = {str(cutoff) + "_shuffle": thisGRN_dict[str(cutoff) + "_shuffle"]}
        
        thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws) 

        my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
        shuffle_stats.append(my_selection)

    combined_df = pd.concat(shuffle_stats, ignore_index = True)

    combined_df.to_csv("output/benchmark_out/FIMO_homoparam5_cutoff" + str(cutoff) + "_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 
    



