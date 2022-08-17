import numpy as np 
import matplotlib.pyplot as plt
from numba import njit
from util import *

'''
SSJRB simulation model. Three parts:
Reservoir releases, gains (into Delta), and Delta pumping
Each component has two functions: step() and fit()
Then all components are combined in the simulate() function
Numba compilation requires only numpy arrays, no pandas objects
'''

@njit 
def reservoir_step(x, dowy, Q, S, K, R_avg, S_avg):
  '''
  Advances reservoir storage from one timestep to the next

    Parameters:
      x (np.array): Reservoir rule parameters (5)
      dowy (int): Day of water year
      Q (float): Inflow, cfs
      S (float): Storage, acre-feet
      K (float): Storage capacity, acre-feet
      R_avg (float): Median release for this day of the year, cfs
      S_avg (float): Median storage for this day of the year, acre-feet

    Returns:
      tuple(float, float): the updated release (cfs) and storage (af)

  '''
  R_target = R_avg * (S / S_avg) ** x[0] # exponential hedging rule
  S_target = S + Q - R_target # assumes 1-day forecast

  excess = S - S_avg * x[4] # volume above top of conservation storage
  if excess > 0 and x[1] < dowy < x[2]: # flood control
    R_target += x[3] * excess
    S_target -= x[3] * excess

  if S_target > K: # spill
    R_target += S_target - K
  elif S_target < 0: # no negative storage
    R_target += S_target

  return (R_target, np.max(np.array([S + Q - R_target, 0])))


@njit
def reservoir_fit(x, dowy, Q, K, Q_avg, R_avg, R_obs, S_avg, S_obs):
  '''
  Evaluate reservoir model against historical observations for a set of parameters 

    Parameters:
      x (np.array): Reservoir rule parameters (5)
      dowy (np.array(int)): Day of water year over the simulation
      Q (np.array(float)): Inflow, cfs
      S (np.array(float)): Storage, acre-feet
      K (float): Storage capacity, acre-feet
      R_avg (np.array(float)): Median release for each day of the year, cfs
      R_obs (np.array(float)): Observed historical release, cfs
      S_avg (np.array(float)): Median storage for each day of the year, acre-feet
      S_obs (np.array(float)): Observed historical storage, acre-feet. 
        Not used currently, but could fit parameters to this instead.

    Returns:
      (float): the negative r**2 value of reservoir releases, to be minimized

  '''
  T = dowy.size
  R,S = np.zeros(T), np.zeros(T)
  S[0] = K / 2 # assume initial condition is half full

  for t in range(1,T):
    inputs = (dowy[t], Q[t], S[t-1], K, R_avg[dowy[t]], S_avg[dowy[t-1]])
    R[t], S[t] = reservoir_step(x, *inputs)
    
  return -np.corrcoef(R_obs, R)[0,1]**2


@njit
def gains_step(x, dowy, Q_total, Q_total_avg, S_total_pct, Gains_avg):
  '''
  Compute gains into the Delta for one timestep

    Parameters:
      x (np.array): Gains parameters (2)
      dowy (int): Day of water year
      Q_total (float): Total inflow to all reservoirs, cfs
      Q_total_avg (float): Average total inflow for this day of the year, cfs
      S_total_pct (float): System-wide reservoir storage, % of median
      Gains_avg (float): Median gains for this day of the year, cfs

    Returns:
      (float): Gains into the Delta for one timestep, cfs

  '''
  G = Gains_avg * S_total_pct ** x[0] # adjust up/down for wet/dry conditions
  if Q_total > Q_total_avg: # high inflows correlated with other tributaries
    G += Q_total * x[1]
  return G


@njit
def gains_fit(x, dowy, Q_total, Q_total_avg, S_total_pct, Gains_avg, Gains_obs):
  '''
  Evaluate Delta gains model against historical observations for a set of parameters

    Parameters:
      x (np.array): Gains parameters (2)
      dowy (np.array(int)): Day of water year
      Q_total (np.array(float)): Total inflow to all reservoirs, cfs
      Q_total_avg (np.array(float)): Average total inflow 
                                     for this day of the year, cfs
      S_total_pct (np.array(float)): System-wide reservoir storage, % of median
      Gains_avg (np.array(float)): Median gains for each day of the year, cfs
      Gains_obs (np.array(float)): Observed historical gains, cfs

    Returns:
      (float): the negative r**2 value of Delta gains, to be minimized

  '''
  T = dowy.size
  G = np.zeros(T)

  for t in range(T):
    inputs = (dowy[t], Q_total[t], Q_total_avg[dowy[t]], S_total_pct[t], Gains_avg[dowy[t]])
    G[t] = gains_step(x, *inputs)

  return -np.corrcoef(Gains_obs, G)[0,1]**2


@njit
def delta_step(x, dowy, Q_in, Pump_pct_avg, S_total_pct):
  '''
  Compute total Delta pumping (Banks + Tracy) for one timestep

    Parameters:
      x (np.array): Delta pumping parameter (1)
      dowy (int): Day of water year
      Q_in (float): Total inflow to the Delta, cfs 
                    (sum of all reservoir outflows plus gains)
      Pump_pct_avg (float): (Average pumping / Average inflow) 
                    for this day of the year, unitless
      S_total_pct (float): System-wide reservoir storage, % of median

    Returns:
      tuple(float, float): Pumping (cfs) and outflow (cfs)

  '''
  P = np.max(np.array([Q_in, 0.0])) * Pump_pct_avg * (S_total_pct ** x[0])
  P = np.min(np.array([P, 11500.0]))
  return (P, Q_in-P)


@njit
def delta_fit(x, dowy, Q_in, Pump_pct_avg, S_total_pct, Pump_obs):
  '''
  Evaluate Delta pumping model against historical observations for a set of parameters

    Parameters:
      x (np.array): Delta pumping parameter (1)
      dowy (np.array(int)): Day of water year
      Q_in (np.array(float)): Total inflow to the Delta, cfs 
                              (sum of all reservoir outflows plus gains)
      Pump_pct_avg (np.array(float)): (Average pumping / Average inflow) 
                                      for each day of the year, unitless
      S_total_pct (np.array(float)): System-wide reservoir storage, % of median
      Pump_obs (np.array(float)): Observed historical pumping, cfs

    Returns:
      (float): the negative r**2 value of Delta pumping, to be minimized

  '''
  T = dowy.size
  P,O = np.zeros(T), np.zeros(T)

  for t in range(T):
    inputs = (dowy[t], Q_in[t], Pump_pct_avg[dowy[t]], S_total_pct[t])
    P[t], O[t] = delta_step(x, *inputs)

  return -np.corrcoef(Pump_obs, P)[0,1]**2


@njit
def simulate(params, K, dowy, Q, Q_avg, R_avg, S_avg, Gains_avg, Pump_pct_avg, DM):
  '''
  Run full system simulation over a given time period.

    Parameters:
      params (tuple(np.array)): Parameter arrays for all reservoirs, gains, and Delta
      K (np.array(float)): Reservoir capacities, acre-feet
      dowy (np.array(int)): Day of water year
      Q (np.array(float, float)): Matrix of inflows at all reservoirs, cfs
      Q_avg (np.array(float, float)): Matrix of median inflows for each reservoir
                                      for each day of the year, cfs
      R_avg (np.array(float, float)): Matrix of median releases for each reservoir
                                      for each day of the year, cfs
      S_avg (np.array(float, float)): Matrix of median storage for each reservoir
                                      for each day of the year, cfs
      Gains_avg (np.array(float)): Median gains for each day of the year, cfs
      Pump_pct_avg (np.array(float)): (Median pumping / median inflow) 
                                      for each day of the year, unitless
      DM (np.array(float)): Demand multiplier, system-wide, unitless
                            (=1.0 for the historical scenario)

    Returns:
      (tuple(np.array, np.array, np.array)): Matrices of timeseries results (cfs) 
        for reservoir releases, storage, and Delta Gains/Inflow/Pumping/Outflow

  '''  
  T,NR = Q.shape
  R,S,G,I,P,O = (np.zeros((T,NR)), np.zeros((T,NR)), 
                 np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T))
  Q_total = Q.sum(axis=1) / cfs_to_afd
  Q_total_avg = Q_avg.sum(axis=1) / cfs_to_afd
  S_total_avg = S_avg.sum(axis=1)
  S[0,:] = K/2

  for t in range(1,T):
    d = dowy[t]

    # 1. Reservoir policies
    for r in range(NR):
      res_demand = R_avg[d,r] * DM[t] # median historical release * demand multiplier
      inputs_r = (d, Q[t,r], S[t-1,r], K[r], res_demand, S_avg[d,r])
      R[t,r], S[t,r] = reservoir_step(params[r], *inputs_r)

    # 2. Gains into Delta
    S_total_pct = S[t].sum() / S_total_avg[d]
    inputs_g = (d, Q_total[t], Q_total_avg[d], 
              S_total_pct, np.min(np.array([Gains_avg[d] * DM[t], Gains_avg[d]])))
    G[t] = gains_step(params[NR], *inputs_g)

    # 3. Delta pumping policy
    I[t] = R[t].sum() / cfs_to_afd + G[t]
    inputs_p = (dowy[t], I[t], Pump_pct_avg[d], S_total_pct)
    P[t], O[t] = delta_step(params[NR+1], *inputs_p)

  Delta = np.vstack((G,I,P,O)).T # Gains, Inflow, Pumping, Outflow
  return (R, S, Delta)