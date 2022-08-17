import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
import model
from scipy.optimize import differential_evolution as DE
from util import *
import time 

nodes = json.load(open('data/nodes.json'))
rk = [k for k in nodes.keys() if nodes[k]['type'] == 'reservoir']
K = np.array([nodes[k]['capacity_taf'] * 1000 for k in rk])
medians = pd.read_csv('data/historical_medians.csv', index_col=0)
params = json.load(open('data/params.json'))
params = tuple(np.array(v) for k,v in params.items())

cmip5_scenarios = pd.read_csv('data/cmip5/scenario_names.csv').name.to_list()
lulc_scenarios = pd.read_csv('data/lulc/scenario_names.csv').name.to_list()
df_sim = None

for sc in cmip5_scenarios:
  for sl in lulc_scenarios:

    print(sc,sl)
    df_sim = pd.DataFrame(index=df_Q.index)
    df_Q = pd.read_csv('data/cmip5/%s.csv.zip' % sc, index_col=0, parse_dates=True, compression='zip')
    df_demand = pd.read_csv('data/lulc/%s.csv.zip' % sl, index_col=0, parse_dates=True)

    input_data = get_simulation_data(rk, df_Q, medians, df_demand)
    R,S,Delta = model.simulate(params, K, *input_data)
    
    df_sim['dowy'] = df_Q.dowy

    for i,r in enumerate(rk):
      df_sim[r+'_outflow_cfs'] = R[:,i]
      df_sim[r+'_storage_af'] = S[:,i]

    for i,k in enumerate(['delta_gains_cfs', 'delta_inflow_cfs', 'total_delta_pumping_cfs', 'delta_outflow_cfs']):
      df_sim[k] = Delta[:,i]

    objs = results_to_annual_objectives(df_sim, medians, nodes, rk, df_demand)
    objs.to_csv('output/scenarios/obj_%s_%s.csv.zip' % (sc,sl), compression='zip')
  
    # too big to save  
    # df_sim.to_csv('output/scenarios/sim_%s_%s.csv.zip' % (sc,sl), compression='zip')
