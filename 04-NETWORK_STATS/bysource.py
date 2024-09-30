#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 13:13:45 2024

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

benchRAPToR = pd.read_csv("output/benchmark_DEstats_RAPToR.txt",
                  sep = '\t',
                  index_col = 0)


benchRAW = pd.read_csv("output/benchmark_DEstats_RAW.txt",
                  sep = '\t',
                  index_col = 0)

#%%
obs1 = pd.read_csv("output/benchmark_observations.txt",
                  sep = '\t',
                     index_col=0,
                     encoding='unicode_escape')

#%%
    
mat = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allthree_equalweights.txt")

# combo3_GRN_dict["combo3_FIMOless"] = pd.read_table("output/GRNs/allthree_fimolessweight.txt")
# combo3_GRN_dict["combo3_chipmore"] = pd.read_table("output/GRNs/allthreeChIP2000_chipmoreweight.txt")


# combo3_GRN_dict["FIMO1500ChIPnocut_weighted_unequal"] = pd.read_table("output/GRNs/FIMO1500ChIPnocutexcl_weighted_unequal.txt")

decouple_kws={
    'methods' : ['mlm'],
    'consensus': False,
}
    

#%%

combo3_output_raptor = dc.benchmark(benchRAPToR, obs1, mat, metrics = ['auprc', 'auroc', 'recall'], by = 'source', min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

combo3_output_raw = dc.benchmark(benchRAW, obs1, mat, metrics = ['auprc', 'auroc', 'recall'], by = 'source', min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

#%%

combo3_output_raptor.to_csv("~/Cel_GRN_revisions/output/benchmark_out/recall_bysource.tsv", sep='\t', index=False)  

combo3_output_raw.to_csv("~/Cel_GRN_revisions/output/benchmark_out/recall_bysource_RAW.tsv", sep='\t', index=False)  