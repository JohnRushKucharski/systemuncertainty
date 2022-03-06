# from operator import index
# import numpy as np 
# import matplotlib.pyplot as plt
import pandas as pd
# import json
import os
from pathlib import Path

cfs_to_afd = 2.29568411*10**-5 * 86400

path = f'{Path(__file__).resolve().parent}' # this .py file directory
inpaths = [f'{path}/error_model/simulated.files.outflows', f'{path}/error_model/simulated.files.pumping']
outpaths = [f'{path}/sobol_analysis/ensemble-simulations-flood', f'{path}/sobol_analysis/ensemble-simulations-pumping']

for i in range(len(inpaths)):
  inputfiles = [f for f in os.listdir(inpaths[i]) if os.path.isfile(os.path.join(inpaths[i], f)) and f.endswith('.csv')]
  for fname in inputfiles:
    df = pd.read_csv(f'{inpaths[i]}/{fname}', index_col=0, parse_dates=True).set_index('Date')
    df.index = pd.to_datetime(df.index)
    sname = fname.split('.csv')[0].replace('.','-')
    if i == 0: # flood
      df = df.resample('AS-OCT').max()
      df.to_csv(f'{outpaths[i]}/{sname}_annmax.csv')
    else: # pumping
      df = df.resample('AS-OCT').sum() * cfs_to_afd / 1000 # taf
      df.to_csv(f'{outpaths[i]}/{sname}_annsum.csv')
    

#cmip5_scenarios = pd.read_csv('../data/cmip5/scenario_names.csv').name.to_list()
#cfs_to_afd = 2.29568411*10**-5 * 86400

# first save zips
# for sc in cmip5_scenarios:
#   sc_mod = sc.replace('-','.')
#   print(sc_mod)
#   # f = '/Users/jon/scratch/scott-ensemble-simulation-unzip/%s.csv' % sc_mod
#   f = '/Users/jon/scratch/scott-ensemble-simulation-pumping-unzip/%s.csv' % sc_mod
#   df = pd.read_csv(f, index_col=0, parse_dates=True)

#   # df.to_csv('scott-ensemble-simulations/%s.csv.zip' % sc, compression='zip')
#   df.to_csv('scott-ensemble-simulations-pumping/%s.csv.zip' % sc, compression='zip')

# then save annual versions
# for sc in cmip5_scenarios:
#   print(sc)

#   # flooding - maximum annual cfs
#   f = 'scott-ensemble-simulations-flood/%s.csv.zip' % sc
#   df = pd.read_csv(f, index_col='Date', parse_dates=True)
#   df.drop('Unnamed: 0', axis=1, inplace=True)
#   df = df.resample('AS-OCT').max()
#   df.to_csv('scott-ensemble-simulations-flood/%s_annmax.csv' % sc)

#   # pumping - annual sum TAF
#   f = 'scott-ensemble-simulations-pumping/%s.csv.zip' % sc
#   df = pd.read_csv(f, index_col='Date', parse_dates=True)
#   df.drop('Unnamed: 0', axis=1, inplace=True)
#   df = df.resample('AS-OCT').sum() * cfs_to_afd / 1000
#   df.to_csv('scott-ensemble-simulations-pumping/%s_annsum.csv' % sc)