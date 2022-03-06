import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path

wd = f'{Path(__file__).resolve().parent}'
cfs_to_afd = 2.29568411*10**-5 * 86400

for obj in ['flood', 'pumping']:
  for gcms in ['allgcms', 'real5', 'priority4']:
    Y = np.loadtxt(f'{wd}/sobol_{gcms}_outputs_{obj}.csv', delimiter=',')
    #Y = np.loadtxt('sobol_%s_outputs_%s.csv' % (gcms,obj), delimiter=',')
    t = np.arange(2020,2099)
    df = pd.DataFrame(Y.T, index=pd.date_range('2020','2099',freq='AS-OCT'))
    print(df)

    if obj == 'flood':
      df_hist = pd.read_csv(f'{wd}/obj_historical.csv', parse_dates=True, index_col=0).Delta_Peak_Inflow_cfs
    else:
      df_hist = pd.read_csv(f'{wd}/sim_historical.csv', parse_dates=True, index_col=0).total_delta_pumping_cfs
      df_hist = df_hist.resample('AS-OCT').sum() * cfs_to_afd / 1000

    for a in ['raw', 'avg']:
              
      w = 1 if a == 'raw' else 10

      df = df.rolling(w).mean()
      ax = df.median(axis=1).plot(color='k')
      df.quantile(0.9, axis=1).plot(color='steelblue', ax=ax)
      df.max(axis=1).plot(color='lightsteelblue', ax=ax)
      df.min(axis=1).plot(color='lightsteelblue', ax=ax)
      df.quantile(0.1, axis=1).plot(color='steelblue', ax=ax)
      df_hist.rolling(w).mean().plot(ax=ax, color='r')
      plt.legend(['Median', 'Max', 'Min', '90%', '10%', 'Historical'])
      ax.set_title('%s objective (%s)' % (obj, a))

      plt.savefig(f'{wd}/figs/future_ranges/svg/{obj}_{gcms}_{a}.svg')
      plt.savefig(f'{wd}/figs/future_ranges/png/{obj}_{gcms}_{a}.png', dpi=300)
      #plt.savefig('figs/future_ranges/svg/%s_%s_%s.svg' % (obj,gcms,a))
      #plt.savefig('figs/future_ranges/png/%s_%s_%s.png' % (obj,gcms,a), dpi=300)
      plt.close()