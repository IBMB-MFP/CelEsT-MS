#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 19:42:29 2024

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

cutoffs_vec = [500, 1000, 1500, 2000, 2500, 3000, 5000, 10000]
fdr_vec = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

#%%    

for cutoff in cutoffs_vec:
    
    for fdr in fdr_vec: 
    
        file_name = "~/Cel_GRN_manuscript/output/GRNs/TF_orthprobs_cut" + str(cutoff) + "_fdr" + str(fdr) + ".txt"
        
        thisGRN_dict = {}
        
        thisGRN_dict[str(cutoff) + "_fdr" + str(fdr)] = pd.read_table(file_name)  
            
        thisGRN_kws = {
            
            str(cutoff) + "_fdr" + str(fdr):{
                'methods' : ['mlm', 'ulm', 'wsum'],
                'consensus' : True,
                # 'dense' : True
                }
            }
            
        thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
        
        # dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
        
        thisoutput.to_csv("output/benchmark_out/TF_orthprobs_cut" + str(cutoff) + "_fdr" + str(fdr) + "_benchRAPToR.tsv", sep='\t', index=False)  
      
