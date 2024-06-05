#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:39:15 2024

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
    
SJARACNE_GRN_dict = {}

SJARACNE_GRN_dict["SJARACNE_GRN"] = pd.read_table("output/GRNs/SJARACNEage_GRN.txt")  # Assuming you're using pandas to read the file
    
SJARACNE_GRN_kws = {
    
    "SJARACNE_GRN":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },

    }
    
SJARACNE_output_raw = dc.benchmark(benchRAW, obs1, SJARACNE_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = SJARACNE_GRN_kws)  # Assuming myfunction is defined elsewhere

dc.plot_metrics_scatter_cols(SJARACNE_output_raw, col = 'method', figsize = (9, 5), groupby = 'net')

SJARACNE_output_raw.to_csv("output/benchmark_out/SJARACNE_benchRAW.tsv", sep='\t', index=False)  # Writing to a file

#%%
    
SJARACNE_output_raptor = dc.benchmark(benchRAPToR, obs1, SJARACNE_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = SJARACNE_GRN_kws)  # Assuming myfunction is defined elsewhere

dc.plot_metrics_scatter_cols(SJARACNE_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

SJARACNE_output_raptor.to_csv("output/benchmark_out/SJARACNE_benchRAPToR.tsv", sep='\t', index=False)  # Writing to a file

#%%

raptor_shuffle_stats = []

for i in range(1, 101):

    SJARACNE_GRN_dict = {}

    SJARACNE_GRN_dict["SJARACNE_GRN"] = pd.read_table("output/GRNs/SJARACNEage_GRN.txt")  # Assuming you're using pandas to read the file
    
    SJARACNE_GRN_dict["SJARACNE_GRN_shuffle"] = dc.shuffle_net(net = SJARACNE_GRN_dict["SJARACNE_GRN"], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
    
    shuffledict = {"SJARACNE_GRN_shuffle": SJARACNE_GRN_dict["SJARACNE_GRN_shuffle"]}
    
    thisoutput = dc.benchmark(benchRAPToR, obs1, shuffledict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = SJARACNE_GRN_kws)  # Assuming myfunction is defined elsewhere

    my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
    
    raptor_shuffle_stats.append(my_selection)

combined_df_raptor = pd.concat(raptor_shuffle_stats, ignore_index = True)

combined_df_raptor.to_csv("output/benchmark_out/SJARACNE_GRN_benchraptor_shufflestats.tsv", sep = '\t', index=False)  # Set index=False to exclude the index from the output


