#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:41:47 2024

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

family_vec = ['ZF - NHR',
              'ZF - C2H2',
              'HD',
              'bHLH',
              'WH',
              'bZIP']

for family in family_vec:

    tempfambench = pd.read_csv("output/" + str(family) + "_family_benchRAPToR.txt", sep = '\t', index_col = 0)
        
    tempfamobs = pd.read_csv("output/" + str(family) + "_family_obs.txt", sep = '\t', index_col = 0, encoding='unicode_escape')

    thisGRN_dict = {}
        
    thisGRN_dict['CelEsT'] = pd.read_table("output/GRNs/allthree_equalweights.txt")  
            
    thisGRN_kws = {
            
        'CelEsT':{
            'methods' : ['mlm'],
            'consensus' : True
            # 'dense' : True
            }
        }
        
    thisoutput = dc.benchmark(tempfambench, tempfamobs, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
        
    # dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
        
    thisoutput.to_csv("output/benchmark_out/FAMILY_" + str(family) + "_benchRAPToR.tsv", sep='\t', index=False)  
      

#%%

for family in family_vec:

    tempfambench = pd.read_csv("output/" + str(family) + "_family_benchRAW.txt", sep = '\t', index_col = 0)
        
    tempfamobs = pd.read_csv("output/" + str(family) + "_family_obs.txt", sep = '\t', index_col = 0, encoding='unicode_escape')

    thisGRN_dict = {}
        
    thisGRN_dict['CelEsT'] = pd.read_table("output/GRNs/allthree_equalweights.txt")  
            
    thisGRN_kws = {
            
        'CelEsT':{
            'methods' : ['mlm'],
            'consensus' : True
            # 'dense' : True
            }
        }
        
    thisoutput = dc.benchmark(tempfambench, tempfamobs, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  
        
    # dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
        
    thisoutput.to_csv("output/benchmark_out/FAMILY_" + str(family) + "_benchRAW.tsv", sep='\t', index=False)  
   
#%%

thisGRN_kws = {
            
        'CelEsT':{
            'methods' : ['mlm'],
            'consensus' : True
            # 'dense' : True
            },
        'CelEsT_shuffle':{
            'methods' : ['mlm'],
            'consensus' : True
            }
        }


thisGRN_dict = {}
        
thisGRN_dict['CelEsT'] = pd.read_table("output/GRNs/allthree_equalweights.txt")  

for family in family_vec:
     
    family_raw_shuffle_stats = []
    
    tempfambench = pd.read_csv("output/" + str(family) + "_family_benchRAW.txt", sep = '\t', index_col = 0)
        
    tempfamobs = pd.read_csv("output/" + str(family) + "_family_obs.txt", sep = '\t', index_col = 0, encoding='unicode_escape')
          
    for i in range(1, 101):
        
        thisGRN_dict['CelEsT_shuffle'] = dc.shuffle_net(net = thisGRN_dict['CelEsT'], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
      
        thisoutput = dc.benchmark(tempfambench, tempfamobs, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = thisGRN_kws)

        my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
                
        family_raw_shuffle_stats.append(my_selection)

    combined_df_raptor = pd.concat(family_raw_shuffle_stats, ignore_index = True)

    combined_df_raptor.to_csv("output/benchmark_out/" + str(family) + "_benchRAW_shufflestats.tsv", sep = '\t', index=False) 
 

#%%

thisGRN_kws = {
            
        'CelEsT':{
            'methods' : ['mlm'],
            'consensus' : True
            # 'dense' : True
            },
        'CelEsT_shuffle':{
            'methods' : ['mlm'],
            'consensus' : True
            }
        }


thisGRN_dict = {}
        
thisGRN_dict['CelEsT'] = pd.read_table("output/GRNs/allthree_equalweights.txt")  

for family in family_vec:
     
    family_raw_shuffle_stats = []
    
    tempfambench = pd.read_csv("output/" + str(family) + "_family_benchRAPToR.txt", sep = '\t', index_col = 0)
        
    tempfamobs = pd.read_csv("output/" + str(family) + "_family_obs.txt", sep = '\t', index_col = 0, encoding='unicode_escape')
          
    for i in range(1, 101):
        
        thisGRN_dict['CelEsT_shuffle'] = dc.shuffle_net(net = thisGRN_dict['CelEsT'], target='target', weight='weight', seed = i).drop_duplicates(['source', 'target'])
      
        thisoutput = dc.benchmark(tempfambench, tempfamobs, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws = thisGRN_kws)

        my_selection = thisoutput[thisoutput['metric'].isin(["auroc", "auprc"])]
                
        family_raw_shuffle_stats.append(my_selection)

    combined_df_raptor = pd.concat(family_raw_shuffle_stats, ignore_index = True)

    combined_df_raptor.to_csv("output/benchmark_out/" + str(family) + "_benchRAPToR_shufflestats.tsv", sep = '\t', index=False) 
 