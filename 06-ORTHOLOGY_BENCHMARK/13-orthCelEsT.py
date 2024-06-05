#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 09:51:49 2024

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

orthCelEsT = pd.read_csv("output/GRNs/orthCelEsT_equalweights.txt",
                  sep = '\t')

CelEsT = pd.read_csv("output/GRNs/allthree_equalweights.txt",
                     sep = "\t")

#%%
FIMOnets = {

    'orthCelEsT': orthCelEsT ,
    'CelEsT': CelEsT
   
   }

decouple_kws = {

    'orthCelEsT':{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True
        },
    'CelEsT':{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True
        },
    
    }
    

#%%
othCelEsT_df = dc.benchmark(benchRAPToR, obs1, FIMOnets, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = decouple_kws)

dc.plot_metrics_scatter_cols(othCelEsT_df, col = 'method', figsize = (9, 5), groupby = 'net')

othCelEsT_df.to_csv("output/benchmark_out/orthCelEsT_df.tsv", sep='\t', index=False)  

#%%

for key in FIMOnets:
     
    othCelEsT_raptor_shuffle_stats = []

    for i in range(1, 101):

        temp_GRN_dict = {}

        temp_GRN_dict[str(key)] = FIMOnets[str(key)]
        
        temp_GRN_dict["temp_GRN_shuffle"] = dc.shuffle_net(net = temp_GRN_dict[str(key)], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
        
        shuffledict = {"temp_GRN_shuffle": temp_GRN_dict["temp_GRN_shuffle"]}
        
        thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = decouple_kws)

        my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
        
        othCelEsT_raptor_shuffle_stats.append(my_selection)

    TForth_shufflestat_raptor = pd.concat(othCelEsT_raptor_shuffle_stats, ignore_index = True)

    TForth_shufflestat_raptor.to_csv("output/benchmark_out/" + str(key) + "_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 
 
    
    
