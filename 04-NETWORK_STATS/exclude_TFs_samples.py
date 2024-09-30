#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:38:23 2024

@author: mfperez
"""


#%%
import pandas as pd
import decoupler as dc
import numpy as np
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

decouple_kws={
    'methods' : ['mlm'],
    'consensus': False,
}

obs1_limit = obs1[obs1['ExperimentNo'].isin(benchRAPToR.index)]

uniqueTFs = obs1_limit['target_gseq'].unique()

for TF in uniqueTFs:

    obs1_tempexcl = obs1_limit[obs1_limit['target_gseq'] != TF]
    
    benchRAPToR_tempexcl = benchRAPToR[benchRAPToR.index.isin(obs1_tempexcl['ExperimentNo'])]

    temp_benchout = dc.benchmark(benchRAPToR_tempexcl, obs1_tempexcl, mat, metrics = ['auprc', 'auroc'], min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

    temp_benchout.to_csv("~/Cel_GRN_revisions/output/benchmark_out/EXCLUDE_" + str(TF) + ".tsv", sep='\t', index=False)  

#%%

mat = pd.read_table("~/Cel_GRN_manuscript/output/GRNs/allthree_equalweights.txt")

decouple_kws={
    'methods' : ['mlm'],
    'consensus': False,
}

obs1_limit = obs1[obs1['ExperimentNo'].isin(benchRAPToR.index)]

uniqueTFs = obs1_limit['target_gseq'].unique()

for TF in uniqueTFs:

    obs1_tempexcl = obs1_limit[obs1_limit['target_gseq'] != TF]
    
    benchRAPToR_tempexcl = benchRAPToR[benchRAPToR.index.isin(obs1_tempexcl['ExperimentNo'])]

    temp_benchout = dc.benchmark(benchRAPToR_tempexcl, obs1_tempexcl, mat, metrics = ['auprc', 'auroc'], min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

    temp_benchout.to_csv("~/Cel_GRN_revisions/output/benchmark_out/EXCLUDE_" + str(TF) + ".tsv", sep='\t', index=False)  


#%%

mat = pd.read_table("output/GRNs/allthree_equalweights.txt")

decouple_kws={
    'methods' : ['mlm'],
    'consensus': False,
}

obs1_limit = obs1[obs1['ExperimentNo'].isin(benchRAW.index)]
obs1_limit = obs1_limit[obs1_limit['target_gseq'].isin(mat['source'])]

for i in range(1, 1001):

    TFsample = np.random.choice(obs1_limit['target_gseq'].unique(), size=10, replace=False)

    obs1_tempexcl = obs1_limit[~(obs1_limit['target_gseq'].isin(TFsample))]
    
    benchRAW_tempexcl = benchRAW[benchRAW.index.isin(obs1_tempexcl['ExperimentNo'])]

    temp_benchout = dc.benchmark(benchRAW_tempexcl, obs1_tempexcl, mat, metrics = ['auprc', 'auroc'], min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

    temp_benchout.to_csv("~/Cel_GRN_revisions/output/benchmark_out/EXCLUDEsample_" + str(i) + "_RAW.tsv", sep='\t', index=False)  

#%%

mat = pd.read_table("output/GRNs/allthree_equalweights.txt")

decouple_kws={
    'methods' : ['mlm'],
    'consensus': False,
}

obs1_limit = obs1[obs1['ExperimentNo'].isin(benchRAW.index)]
obs1_limit = obs1_limit[obs1_limit['target_gseq'].isin(mat['source'])]

shuffle_stats = []

for i in range(1, 1001):

    TFsample = np.random.choice(obs1_limit['target_gseq'].unique(), size=10, replace=False)

    shufflednet = dc.shuffle_net(net = mat, target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
    
    obs1_tempexcl = obs1_limit[~(obs1_limit['target_gseq'].isin(TFsample))]
 
    benchRAW_tempexcl = benchRAW[benchRAW.index.isin(obs1_tempexcl['ExperimentNo'])]

    temp_benchout = dc.benchmark(benchRAW_tempexcl, obs1_tempexcl, shufflednet, metrics = ['auprc', 'auroc'], min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

    my_selection = temp_benchout[temp_benchout['metric'].isin(["auroc", "auprc"])]
    
    shuffle_stats.append(my_selection)

combined_df = pd.concat(shuffle_stats, ignore_index = True)

combined_df.to_csv("output/benchmark_out/EXCLUDEsample_benchRAW_shufflestats.tsv", sep = '\t', index=False)  

#%%

shuffle_stats = []

for i in range(1, 1001):

    TFsample = np.random.choice(obs1_limit['target_gseq'].unique(), size=10, replace=False)

    shufflednet = dc.shuffle_net(net = mat, target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
    
    obs1_tempexcl = obs1_limit[~(obs1_limit['target_gseq'].isin(TFsample))]
 
    benchRAPToR_tempexcl = benchRAPToR[benchRAPToR.index.isin(obs1_tempexcl['ExperimentNo'])]

    temp_benchout = dc.benchmark(benchRAPToR_tempexcl, obs1_tempexcl, shufflednet, metrics = ['auprc', 'auroc'], min_exp = 1, perturb = 'target_gseq', sign = -1, verbose = True, decouple_kws = decouple_kws)  

    my_selection = temp_benchout[temp_benchout['metric'].isin(["auroc", "auprc"])]
    
    shuffle_stats.append(my_selection)

combined_df = pd.concat(shuffle_stats, ignore_index = True)

combined_df.to_csv("output/benchmark_out/EXCLUDEsample_benchRAPToR_shufflestats.tsv", sep = '\t', index=False)  
