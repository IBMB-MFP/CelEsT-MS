#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 20:27:27 2024

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
    
combo_GRN_dict = {}

combo_GRN_dict["combo3"] = pd.read_table("output/GRNs/allthree_equalweights.txt")  
combo_GRN_dict["combo3negreg"] = pd.read_table("output/GRNs/allthree_equalweights_negreg.txt")

combo_GRN_kws = {
    
    "combo3":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },
    "combo3negreg":{
        'methods' : ['mlm', 'ulm', 'wsum'],
        'consensus' : True,
        # 'dense' : True
        },

    }

#%%

combo_output_raptor = dc.benchmark(benchRAPToR, obs1, combo_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo_GRN_kws)  

dc.plot_metrics_scatter_cols(combo_output_raptor, col = 'method', figsize = (9, 5), groupby = 'net')

combo_output_raptor.to_csv("output/benchmark_out/negregtest.txt")

#%%

combo_output_raw = dc.benchmark(benchRAW, obs1, combo_GRN_dict, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = combo_GRN_kws)  

dc.plot_metrics_scatter_cols(combo_output_raw, col = 'method', figsize = (9, 5), groupby = 'net')

