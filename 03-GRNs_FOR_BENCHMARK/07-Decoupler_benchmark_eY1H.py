#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:40:40 2024

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
    
e1YH_GRN_dict = {}

e1YH_GRN_dict["e1YH_GRN"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/walhout_highqual_cutoff15.txt")  
    
e1YH_GRN_kws = {
    
    "e1YH_GRN":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },

    }
    
# e1YH_output_raw = dc.benchmark(benchRAW, obs1, e1YH_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = e1YH_GRN_kws)  

# dc.plot_metrics_scatter_cols(e1YH_output_raw, col = 'method', figsize = (9, 5), groupby = 'net')

# e1YH_output_raw.to_csv("output/benchmark_out/e1YH_benchRAW.tsv", sep='\t', index=False)  

#%%

e1YH_output_raptor = dc.benchmark(benchRAPToR, obs1, e1YH_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = e1YH_GRN_kws)  

dc.plot_metrics_scatter_cols(e1YH_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

e1YH_output_raptor.to_csv("output/benchmark_out/e1YH_benchRAPToR.tsv", sep='\t', index=False)  

#%%

raptor_shuffle_stats = []

for i in range(1, 101):

    e1YH_GRN_dict = {}

    e1YH_GRN_dict["e1YH_GRN"] = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/walhout_highqual_cutoff15.txt")
    
    e1YH_GRN_dict["e1YH_GRN_shuffle"] = dc.shuffle_net(net = e1YH_GRN_dict["e1YH_GRN"], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
    
    shuffledict = {"e1YH_GRN_shuffle": e1YH_GRN_dict["e1YH_GRN_shuffle"]}
    
    thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = e1YH_GRN_kws)

    my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
    
    raptor_shuffle_stats.append(my_selection)

combined_df_raptor = pd.concat(raptor_shuffle_stats, ignore_index = True)

combined_df_raptor.to_csv("output/benchmark_out/e1YH_GRN_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 

#%%

# raw_shuffle_stats = []

# for i in range(1, 101):

#     e1YH_GRN_dict = {}

#     e1YH_GRN_dict["e1YH_GRN"] = pd.read_table("output/GRNs/walhout_highqual_cutoff15.txt") 
    
#     e1YH_GRN_dict["e1YH_GRN_shuffle"] = dc.shuffle_net(net = e1YH_GRN_dict["e1YH_GRN"], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
    
#     shuffledict = {"e1YH_GRN_shuffle": e1YH_GRN_dict["e1YH_GRN_shuffle"]}
    
#     thisoutput = dc.benchmark(benchRAW, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = e1YH_GRN_kws) 

#     my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
    
#     raw_shuffle_stats.append(my_selection)

# combined_df_raw = pd.concat(raw_shuffle_stats, ignore_index = True)

# combined_df_raw.to_csv("output/benchmark_out/e1YH_GRN_benchRAW_shufflestats.tsv", sep = '\t', index=False) 
