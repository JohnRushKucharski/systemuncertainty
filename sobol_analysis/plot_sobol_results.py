import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path

wd = f'{Path(__file__).resolve().parent}'

for obj in ['flood', 'pumping']:
  for gcms in ['allgcms', 'real5', 'priority4']:

    df = pd.read_csv(f'{wd}/sobol_{gcms}_si_{obj}.csv', index_col=0, parse_dates=True)
    #df = pd.read_csv('sobol_%s_si_%s.csv' % (gcms,obj), index_col=0, parse_dates=True)
    cond = df.columns.str.contains('_T')

    for s in ['si', 'st']:
      df_s = df.loc[:, ~cond] if s == 'si' else df.loc[:, cond]

      for a in ['raw', 'avg']:
        
        if a == 'raw':
          df_s.plot.area()
        else:
          df_s.rolling(10).mean().plot.area()
          # plt.xlim(['2030','2100'])

        plt.savefig(f'{wd}/figs/sobol/svg/{obj}_{gcms}_{s}_{a}.svg')
        plt.savefig(f'{wd}/figs/sobol/png/{obj}_{gcms}_{s}_{a}.png', dpi=300)
        # plt.savefig('figs/sobol/svg/%s_%s_%s_%s.svg' % (obj,gcms,s,a))
        # plt.savefig('figs/sobol/png/%s_%s_%s_%s.png' % (obj,gcms,s,a), dpi=300)
        plt.close()