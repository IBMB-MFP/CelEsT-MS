#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 14:24:18 2024

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

cutoffs_vec = [100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000]

for cutoff in cutoffs_vec:
    
    file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTexcl.txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  # Assuming you're using pandas to read the file
    
    thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
    
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        }
        
    thisoutput = dc.benchmark(benchRAW, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  # Assuming myfunction is defined elsewhere
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTexcl_benchRAW_dense.tsv", sep='\t', index=False)  # Writing to a file

#%%    

for cutoff in cutoffs_vec:
    
    file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTexcl.txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  # Assuming you're using pandas to read the file
    
    thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
    
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        }
        
    thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  # Assuming myfunction is defined elsewhere
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTexcl_benchRAPToR_dense.tsv", sep='\t', index=False)  # Writing to a file
  
#%%  
    
for cutoff in cutoffs_vec:
    
    file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTincl.txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  # Assuming you're using pandas to read the file
    
    thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
    
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        }
        
    thisoutput = dc.benchmark(benchRAW, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  # Assuming myfunction is defined elsewhere
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTincl_benchRAW_dense.tsv", sep='\t', index=False)  # Writing to a file
  
#%%      

for cutoff in cutoffs_vec:
    
    file_name = "output/GRNs/allChIP_" + str(cutoff) + "_HOTincl.txt"
    
    thisGRN_dict = {}
    
    thisGRN_dict[str(cutoff)] = pd.read_table(file_name)  # Assuming you're using pandas to read the file
    
    thisGRN_dict[str(cutoff) + "_shuffle"] = dc.shuffle_net(net = thisGRN_dict[str(cutoff)], target='target', weight='weight').drop_duplicates(['source', 'target'])
    
    thisGRN_kws = {
        
        str(cutoff):{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        str(cutoff) + "_shuffle":{
            'methods' : ['mlm', 'ulm', 'wsum'],
            'consensus' : True,
            'dense' : True
            },
        }
        
    thisoutput = dc.benchmark(benchRAPToR, obs1, thisGRN_dict, perturb='target_gseq', sign=-1, verbose=True, decouple_kws=thisGRN_kws)  # Assuming myfunction is defined elsewhere
    
    dc.plot_metrics_scatter_cols(thisoutput, col = 'method', figsize = (9, 5), groupby = 'net')
    
    thisoutput.to_csv("output/benchmark_out/allChIP_cutoff" + str(cutoff) + "_HOTincl_benchRAPToR_dense.tsv", sep='\t', index=False)  # Writing to a file
    