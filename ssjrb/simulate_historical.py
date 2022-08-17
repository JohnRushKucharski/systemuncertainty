import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
import model
from scipy.optimize import differential_evolution as DE
from util import *
import time 

nodes = json.load(open('data/nodes.json'))
df_hist = pd.read_csv('data/historical.csv', index_col=0, parse_dates=True)['10-01-1997':]
medians = pd.read_csv('data/historical_medians.csv', index_col=0)
params = json.load(open('data/params.json'))
params = tuple(np.array(v) for k,v in params.items())

rk = [k for k in nodes.keys() if nodes[k]['type'] == 'reservoir']
K = np.array([nodes[k]['capacity_taf'] * 1000 for k in rk])

# run simulation and put results in dataframe
# desktop runtime 0.014s with numba, 0.19s without
input_data = get_simulation_data(rk, df_hist, medians)
R,S,Delta = model.simulate(params, K, *input_data)

df_sim = pd.DataFrame(index=df_hist.index)
df_sim['dowy'] = df_hist.dowy

for i,r in enumerate(rk):
  df_sim[r+'_outflow_cfs'] = R[:,i]
  df_sim[r+'_storage_af'] = S[:,i]

for i,k in enumerate(['delta_gains_cfs', 'delta_inflow_cfs', 'total_delta_pumping_cfs', 'delta_outflow_cfs']):
  df_sim[k] = Delta[:,i]

objs = results_to_annual_objectives(df_sim, medians, nodes, rk)
df_sim.to_csv('output/sim_historical.csv')
objs.to_csv('output/obj_historical.csv')