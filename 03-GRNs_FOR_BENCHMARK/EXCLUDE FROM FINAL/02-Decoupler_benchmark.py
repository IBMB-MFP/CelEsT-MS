#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 11:40:42 2024

@author: mfperez
"""

#%%
import pandas as pd
import decoupler as dc
import os

#%%
benchRAW = pd.read_csv("~/Cel_GRN_manuscript/output/benchmark_DEstats_RAW.txt",
                  sep = '\t',
                  index_col = 0)
#%%
obs1 = pd.read_csv("~/Cel_GRN_manuscript/output/benchmark_observations.txt",
                  sep = '\t',
                     index_col=0,
                     encoding='unicode_escape')

#%%

ChIP2000 = pd.read_csv("~/Cel_GRN_manuscript/output/GRNs/allChIP_2000_HOTincl.txt",
                               sep = '\t')

ChIP2000_noHOT = pd.read_csv("~/Cel_GRN_manuscript/output/GRNs/allChIP_2000_HOTexcl.txt",
                               sep = '\t')

ChIP2500 = pd.read_csv("~/Cel_GRN_manuscript/output/GRNs/allChIP_2500_HOTincl.txt",
                               sep = '\t')

ChIP2500_noHOT = pd.read_csv("~/Cel_GRN_manuscript/output/GRNs/allChIP_2500_HOTexcl.txt",
                               sep = '\t')

#%%

# build dictionary of networks to test

test_nets = {

   'ChIP2000': ChIP2000,
   'ChIP2000_noHOT': ChIP2000_noHOT,
   'ChIP2500': ChIP2500,
   'ChIP2500_noHOT': ChIP2500_noHOT,

    }

#%%

decouple_kws = {
    
    'ChIP2000':{
        'methods' : ['mlm', 'ulm'],
        'consensus' : True,
        'dense' : True
        },
    'ChIP2000_noHOT':{
        'methods' : ['mlm', 'ulm'],
        'consensus' : True,
        'dense' : True
        },
    'ChIP2500':{
        'methods' : ['mlm', 'ulm'],
        'consensus' : True,
        'dense' : True
        },
    'ChIP2500_noHOT':{
        'methods' : ['mlm', 'ulm'],
        'consensus' : True,
        'dense' : True
        },
    }


#%%

test_df = dc.benchmark(benchRAW, obs1, test_nets, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=decouple_kws)

dc.plot_metrics_scatter_cols(test_df, col = 'method', figsize = (9, 5), groupby = 'net')


