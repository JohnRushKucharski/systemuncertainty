import gc
import numpy as np
from matplotlib import pyplot as plt
from SALib.sample import saltelli
from SALib.analyze import sobol
import pandas as pd
from pathlib import Path
import os

wd = f'{Path(__file__).resolve().parent}'

# choice of GCMs
# https://www.nature.com/articles/s41598-019-46169-w#MOESM1
# Real-5: access1-0, access1-3 (don't have in USBR), canesm2, cnrm-cm5, gfdl-cm3
# https://www.energy.ca.gov/sites/default/files/2019-11/Projections_CCCA4-CEC-2018-006_ADA.pdf
# Priority 4 CEC: hadgem2-es, canesm2, cnrm-cm5, miroc5

gcm_list = ['access1-0', 'bcc-csm1-1', 'bcc-csm1-1-m', 'canesm2', 'ccsm4', 'cesm1-bgc', 'cesm1-cam5',
            'cmcc-cm', 'cnrm-cm5', 'csiro-mk3-6-0', 'fgoals-g2', 'fio-esm', 'gfdl-cm3', 'gfdl-esm2g',
            'gfdl-esm2m', 'giss-e2-r', 'hadgem2-ao', 'hadgem2-cc',
            'hadgem2-es', 'inmcm4', 'ipsl-cm5a-mr', 'ipsl-cm5b-lr', 'miroc5', 'miroc-esm', 'miroc-esm-chem',
            'mpi-esm-lr', 'mpi-esm-mr', 'mri-cgcm3', 'noresm1-m']

# settings ~~~~~~~~~~~~~~~~~~~~~~
rcp_list = ['rcp45', 'rcp85']
obj = 'flood'
gcms = 'priority4' # allgcms, real5, priority4
rerun = False

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if gcms == 'real5':
  gcm_list = ['access1-0', 'canesm2', 'cnrm-cm5', 'gfdl-cm3']
elif gcms == 'priority4':
  gcm_list = ['hadgem2-es', 'canesm2', 'cnrm-cm5', 'miroc5']

# gcm, rcp, system
problem = {
        'num_vars': 3,
        'names': ['rcp', 'gcm', 'system'],
        'bounds': [[0, 2], [0, 4], [1, 101]]}

# Generate samples
X = saltelli.sample(problem, 8192) # 16384 used for allgcms)
X = X.astype(int)  # scenario numbers are integers
N = len(X)
t = np.arange(2020,2099)

if rerun:
  dir_error_model = f'{Path(__file__).resolve().parent.parent}/error_model/'
  # initiate list of objective values (eventually converted to numpy vector)
  Y = np.zeros((N, len(t)))

  for i,x in enumerate(X):

    if i % 100 == 0: print(i)
    rcp,gcm,sys = rcp_list[x[0]], gcm_list[x[1]], x[2]

    if obj=='flood':
      # does it need the _annmax ?
      fname = f'{dir_error_model}simulated.files.outflows/{gcm}_{rcp}_r1i1p1.csv'
      df = pd.read_csv(fname, index_col=0, parse_dates=True)['X%s' % sys]
    else:
      # does it need the _annsum ?
      fname = f'{dir_error_model}simulated.files.pumping/{gcm}_{rcp}_r1i1p1.csv'
      df = pd.read_csv(fname, index_col=0, parse_dates=True)['X%s' % sys]

    Y[i,:] = df['2020':'2099'].values

  np.savetxt(f'{wd}/sobol_{gcms}_outputs_{obj}.csv', Y, delimiter=',')
  #np.savetxt('sobol_%s_outputs_%s.csv' % (gcms,obj), Y, delimiter=',')

else:
  Y = np.loadtxt(f'{wd}/sobol_{gcms}_outputs_{obj}.csv', delimiter=',')
  #Y = np.loadtxt('sobol_%s_outputs_%s.csv' % (gcms,obj), delimiter=',')

  fname = f'{wd}/sobol_{gcms}_si_{obj}.csv'
  #fname = 'sobol_%s_si_%s.csv' % (gcms,obj)

  # if not os.path.exists(fname):
  df = pd.DataFrame(index=t, columns=['rcp', 'gcm', 'system', 
                                      'rcp:gcm', 'rcp:system', 'gcm:system',
                                      'rcp_T', 'gcm_T', 'system_T'])

  for y in t:
    print(y)
    Si = sobol.analyze(problem, Y[:,y-2020], print_to_console=False)
    df.loc[y, 'rcp'] = Si['S1'][0]
    df.loc[y, 'gcm'] = Si['S1'][1]
    df.loc[y, 'system'] = Si['S1'][2]
    df.loc[y, 'rcp_T'] = Si['ST'][0]
    df.loc[y, 'gcm_T'] = Si['ST'][1]
    df.loc[y, 'system_T'] = Si['ST'][2]
    df.loc[y, 'rcp:gcm'] = Si['S2'][0,1]
    df.loc[y, 'rcp:system'] = Si['S2'][0,2]
    df.loc[y, 'gcm:system'] = Si['S2'][1,2]

  df[df < 0] = 0
  df.to_csv(f'{wd}/sobol_{gcms}_si_{obj}.csv')
  #df.to_csv('sobol_%s_si_%s.csv' % (gcms,obj))

  # else:
  #   df = pd.read_csv(fname, index_col=0, parse_dates=True)

  # df.plot.area()
  # plt.show()
