#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 14:24:18 2024

@author: mfperez
"""


#%%
import pandas as pd
import decoupler as dc
import os

from pathlib import Path

#%%

home_directory = os.path.expanduser("~")

working_directory = os.path.join(home_directory, "Cel_GRN_revisions")

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

# #%%
# for cutoff in cutoffs_vec:
    
#     file_name = "~/Cel_GRN_manuscript/output/GRNs/allChIP_" + str(cutoff) + "_HOTexcl.txt"
    
#     thisGRN_dict = {}
    
#     thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
        
#     thisGRN_kws = {
        
#         str(cutoff):{
#             'methods' : ['mlm', 'ulm', 'wsum'],
#             'consensus' : True,
#             # 'dense' : True
#             }
#         }
        
#     thisoutput = dc.benchmark(benchRAW, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
    
#     dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
#     thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTexcl_benchRAW.tsv", sep='\t', index=False)  

#%%    

for cutoff in cutoffs_vec:
    
    file_name = "~/Cel_GRN_manuscript/output/GRNs/allChIP_" + str(cutoff) + "_HOTexcl.txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
        
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            }
        }
        
    thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTexcl_benchRAPToR.tsv", sep='\t', index=False)  
  
#%%  
    
# for cutoff in cutoffs_vec:
    
#     file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTincl.txt"
    
#     thisGRN_dict = {}
    
#     thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  # Assuming you're using pandas to read the file
        
#     thisGRN_kws = {
        
#         str(cutoff):{
#             'methods' : ['mlm', 'ulm', 'wsum'],
#             'consensus' : True,
#             # 'dense' : True
#             }
#         }
        
#     thisoutput = dc.benchmark(benchRAW, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws) 
    
#     dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
#     thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTincl_benchRAW.tsv", sep='\t', index=False)  
  
#%%      

for cutoff in cutoffs_vec:
    
    file_name = "~/Cel_GRN_manuscript/output/GRNs/allChIP_" + str(cutoff) + "_HOTincl.txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
        
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            # 'dense' : True
            }
        }
        
    thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTincl_benchRAPToR.tsv", sep='\t', index=False)  

#%% RUN FOR HOT REGION CUTOFFS
    
HOTchip_dict = {}

HOTchip_dict["20"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTcut20.txt")  
HOTchip_dict["30"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTcut30.txt")
HOTchip_dict["40"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTcut40.txt")
HOTchip_dict["50"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTcut50.txt")
HOTchip_dict["70"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTcut70.txt")

HOTchip_kws = {
    
    "20":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "30":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "40":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "50":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "70":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        }

    }

#%%

HOTchip_output_raptor = dc.benchmark(benchRAPToR, obs1, HOTchip_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = HOTchip_kws)  

dc.plot_metrics_scatter_cols(HOTchip_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

HOTchip_output_raptor.to_csv("output/benchmark_out/HotCut_bench.tsv", sep='\t', index=False)

#%%     

# run iterations for shuffled networks

for cutoff in cutoffs_vec:

    file_name = "~/Cel_GRN_manuscript/output/GRNs/allChIP_" + str(cutoff) + "_HOTincl.txt"
    
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

    combined_df.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTincl_benchRAPToR_shufflestats.tsv", sep = '\t', index=False)  


# run iterations for shuffled networks

# for cutoff in cutoffs_vec:

#     file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTincl.txt"
    
#     thisGRN_dict = {}
    
#     thisGRN_kws = {
        
#         str(cutoff) + "_shuffle":{
#             'methods' : ['mlm', 'ulm', 'wsum'],
#             'consensus' : True,
#             # 'dense' : True
#             },
#         }
        
#     shuffle_stats = []
    
#     for i in range(1, 101):

#         thisGRN_dict[str(cutoff)] = pd.read_table(file_name) 
    
#         thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
        
#         shuffledict = {str(cutoff) + "_shuffle": thisGRN_dict[str(cutoff) + "_shuffle"]}
        
#         thisoutput = dc.benchmark(benchRAW, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  

#         my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
#         shuffle_stats.append(my_selection)

#     combined_df = pd.concat(shuffle_stats, ignore_index = True)

#     combined_df.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTincl_benchRAW_shufflestats.tsv", sep = '\t', index=False)  

# run iterations for shuffled networks

#%%   

# for cutoff in cutoffs_vec:

#     file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTexcl.txt"
    
#     thisGRN_dict = {}
    
#     thisGRN_kws = {
        
#         str(cutoff) + "_shuffle":{
#             'methods' : ['mlm', 'ulm', 'wsum'],
#             'consensus' : True,
#             # 'dense' : True
#             },
#         }
        
#     shuffle_stats = []
    
#     for i in range(1, 101):

#         thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  
    
#         thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
        
#         shuffledict = {str(cutoff) + "_shuffle": thisGRN_dict[str(cutoff) + "_shuffle"]}
        
#         thisoutput = dc.benchmark(benchRAW, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  

#         my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
#         shuffle_stats.append(my_selection)

#     combined_df = pd.concat(shuffle_stats, ignore_index = True)

#     combined_df.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTexcl_benchRAW_shufflestats.tsv", sep = '\t', index=False)  
    

for cutoff in cutoffs_vec:

    file_name = "~/Cel_GRN_manuscript/output/GRNs/allChIP_" + str(cutoff) + "_HOTexcl.txt"
    
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
        
        thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  # Assuming myfunction is defined elsewhere

        my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
        shuffle_stats.append(my_selection)

    combined_df = pd.concat(shuffle_stats, ignore_index = True)

    combined_df.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTexcl_benchRAPToR_shufflestats.tsv", sep = '\t', index=False)  # Set index=False to exclude the index from the output


