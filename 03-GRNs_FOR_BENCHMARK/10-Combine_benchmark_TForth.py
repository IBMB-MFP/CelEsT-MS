#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:36:43 2024

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
obs2 = pd.read_csv("output/TForthChIPnocutexcl_observations.txt",
                  sep = '\t',
                     index_col=0,
                     encoding='unicode_escape')

#%%
    
combo_GRN_dict = {}

combo_GRN_dict["TForth"] = pd.read_table("output/GRNs/TF_orthprobs_cut1500_fdr0.5.txt")  
combo_GRN_dict["FIMO1500"] = pd.read_table("output/GRNs/FIMO_nohomo_1500.txt")  
combo_GRN_dict["ChIPnocut"] = pd.read_table("output/GRNs/allChIP_10000_HOTexcl.txt")
combo_GRN_dict["TForthChIPnocut_noweights"] = pd.read_table("output/GRNs/TForthChIPnocutexcl_noweights.txt")
combo_GRN_dict["TForthChIPnocut_weighted"] = pd.read_table("output/GRNs/TForthChIPnocutexcl_weighted.txt")
combo_GRN_dict["FIMO1500ChIPnocut_weighted"] = pd.read_table("output/GRNs/FIMO1500ChIPnocutexcl_weighted.txt")
# combo_GRN_dict["TForthChIPnocut_weighted_unequal"] = pd.read_table("output/GRNs/TForthChIPnocutexcl_weighted_unequal.txt")

combo_GRN_kws = {
    
    "TForth":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "FIMO1500":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "ChIPnocut":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "TForthChIPnocut_noweights":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "TForthChIPnocut_weighted":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "FIMO1500ChIPnocut_weighted":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    # "TForthChIPnocut_weighted_unequal":{
    #     'methods' : ['mlm', 'ulm', 'wsum'],
    #     'consensus' : True,
    #     # 'dense' : True
        # },

    }

# #%%

# combo_output_raptor = dc.benchmark(benchRAPToR, obs1, combo_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo_GRN_kws)  

# dc.plot_metrics_scatter_cols(combo_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

#%%

combo_output_raptor = dc.benchmark(benchRAPToR.loc[obs2.index], obs2, combo_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo_GRN_kws)  

dc.plot_metrics_scatter_cols(combo_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

# combo_output_raptor.to_csv("output/benchmark_out/combo_restrictedobs_benchRAPToR.tsv", sep='\t', index=False)

#%%
    
combo3_GRN_dict = {}

combo3_GRN_dict["combo3equal"] = pd.read_table("output/GRNs/allthree_TForth_equalweights.txt")
combo3_GRN_dict["combo3no"] = pd.read_table("output/GRNs/allthree_TForth_noweights.txt")
combo3_GRN_dict["combo3_unequal"] = pd.read_table("output/GRNs/allthree_TForth_weighted_unequal.txt")
combo3_GRN_dict["OLDequal"] = pd.read_table("output/GRNs/allthree_equalweights.txt")

# combo3_GRN_dict["combo3_FIMOless"] = pd.read_table("output/GRNs/allthree_fimolessweight.txt")
# combo3_GRN_dict["combo3_chipmore"] = pd.read_table("output/GRNs/allthreeChIP2000_chipmoreweight.txt")


# combo3_GRN_dict["TForthChIPnocut_weighted_unequal"] = pd.read_table("output/GRNs/TForthChIPnocutexcl_weighted_unequal.txt")

combo3_GRN_kws = {
    
    "combo3equal":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "combo3no0":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "combo3_unequal":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "combo3_OLDequal":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    # "combo3_chipmore":{
    #     'methods' : ['mlm', 'ulm', 'wsum'],
    #     'consensus' : True,
    #     # 'dense' : True
    #     },
    }

# %%

combo3_output_raptor = dc.benchmark(benchRAPToR, obs1, combo3_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo3_GRN_kws)  

dc.plot_metrics_scatter_cols(combo3_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

# %%

combo3_output_raw = dc.benchmark(benchRAW, obs1, combo3_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo3_GRN_kws)  

dc.plot_metrics_scatter_cols(combo3_output_raw, col = 'method', figsize = (9, 5), groupby = 'net')
# combo3_output_raptor.to_csv("output/benchmark_out/combo3_benchRAPToR.tsv", sep='\t', index=False)  

# #%%

# combo3_output_raw = dc.benchmark(benchRAW, obs1, combo3_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo3_GRN_kws)  

# dc.plot_metrics_scatter_cols(combo3_output_raw, col = 'method', figsize = (9, 5), groupby = 'net')

# combo3_output_raw.to_csv("output/benchmark_out/combo3_benchraw.tsv", sep='\t', index=False)  

# #%%

# for key in combo3_GRN_dict:
     
#     combo3_raptor_shuffle_stats = []

#     for i in range(1, 101):

#         temp_GRN_dict = {}

#         temp_GRN_dict[str(key)] = combo3_GRN_dict[str(key)]
        
#         temp_GRN_dict["temp_GRN_shuffle"] = dc.shuffle_net(net = temp_GRN_dict[str(key)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
        
#         shuffledict = {"temp_GRN_shuffle": temp_GRN_dict["temp_GRN_shuffle"]}
        
#         thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = combo3_GRN_kws)

#         my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
#         combo3_raptor_shuffle_stats.append(my_selection)

#     combined_df_raptor = pd.concat(combo3_raptor_shuffle_stats, ignore_index = True)

#     combined_df_raptor.to_csv("output/benchmark_out/" + str(key) + "_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 
 
#     #%%

#     for key in combo3_GRN_dict:
         
#         combo3_raw_shuffle_stats = []

#         for i in range(1, 101):

#             temp_GRN_dict = {}

#             temp_GRN_dict[str(key)] = combo3_GRN_dict[str(key)]
            
#             temp_GRN_dict["temp_GRN_shuffle"] = dc.shuffle_net(net = temp_GRN_dict[str(key)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
            
#             shuffledict = {"temp_GRN_shuffle": temp_GRN_dict["temp_GRN_shuffle"]}
            
#             thisoutput = dc.benchmark(benchRAW, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = combo3_GRN_kws)

#             my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
            
#             combo3_raw_shuffle_stats.append(my_selection)

#         combined_df_raw = pd.concat(combo3_raw_shuffle_stats, ignore_index = True)

#         combined_df_raw.to_csv("output/benchmark_out/" + str(key) + "_benchraw_shufflestats.tsv", sep = '\t', index=False) 
     
#     #%%

#     for key in combo_GRN_dict:
         
#         combo_raptor_shuffle_stats = []

#         for i in range(1, 101):

#             temp_GRN_dict = {}

#             temp_GRN_dict[str(key)] = combo_GRN_dict[str(key)]
            
#             temp_GRN_dict["temp_GRN_shuffle"] = dc.shuffle_net(net = temp_GRN_dict[str(key)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
            
#             shuffledict = {"temp_GRN_shuffle": temp_GRN_dict["temp_GRN_shuffle"]}
            
#             thisoutput = dc.benchmark(benchRAPToR.loc[obs2.index], obs2, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = combo_GRN_kws)

#             my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
            
#             combo_raptor_shuffle_stats.append(my_selection)

#         combined_df_raptor = pd.concat(combo_raptor_shuffle_stats, ignore_index = True)

#         combined_df_raptor.to_csv("output/benchmark_out/" + str(key) + "restrictedobs_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 
     




