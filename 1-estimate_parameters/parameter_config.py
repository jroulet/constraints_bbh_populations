"""Write a configuration file `parameter_grid` with parameter values on 
which `compute_likelihood.py` will generate templates and compute the 
likelihood.
* Usage:
    python parameter_config.py <run_directory/>
"""

import numpy as np
import pandas as pd
import sys
from collections import OrderedDict

try:
    event = sys.argv[1].replace('/', '')
except IndexError:
    sys.exit('Usage:\npython parameter_config.py <run_directory/>')

grid_nums = OrderedDict([('M_chirp', 64 if event != 'GW170608' else 32),
                         ('q'      , 64 if event != 'GW170608' else 128),
                         ('chi_eff', 64 if event != 'GW170608' else 128)
                        ])  # Number of points along each dimension

grid_params = list(grid_nums.keys())
zoom_out = 2.6  # Grid extent relative to the LIGO errorbars

ligo_params = pd.read_csv('{0}/{0}_LIGO_parameters'.format(event), 
                          sep=r'\s+', comment='#')

bounds = {'M_chirp': [1.2, 100],
          'q':       [.07, 1],
          'chi_eff': [-1, 1]
         }

grid_1d, reported = OrderedDict(), OrderedDict()
for par in grid_params:
    reported[par] = ligo_params[par]['Overall']
    min_val = max(reported[par] - zoom_out * ligo_params[par]['Overall_errm'],
                  bounds[par][0])
    max_val = min(reported[par] + zoom_out * ligo_params[par]['Overall_errp'],
                  bounds[par][1])
    grid_1d[par] = np.linspace(min_val, max_val, grid_nums[par])

# Hardcode any that you want to change:
if event == 'GW170608':
    grid_1d['M_chirp'] = np.linspace(8.4, 8.7, grid_nums['M_chirp'])
    grid_1d['chi_eff'] = np.linspace(-.1, .7, grid_nums['chi_eff'])

if event == 'GW170814':
    grid_1d['M_chirp'] = np.linspace(23.5, 29.5, grid_nums['M_chirp'])

if event == 'GW170729':
  grid_1d['M_chirp'] = np.linspace(35, 68, grid_nums['M_chirp'])

if event == 'GW170809':
  grid_1d['M_chirp'] = np.linspace(26, 34, grid_nums['M_chirp'])

if event == 'GW170818':
  grid_1d['M_chirp'] = np.linspace(26, 38, grid_nums['M_chirp'])

if event == 'GW170823':
  grid_1d['M_chirp'] = np.linspace(28, 47, grid_nums['M_chirp'])



# Write grid:
grid = OrderedDict(zip(grid_params, [x.flatten() for x in np.meshgrid(
    *[grid_1d[par] for par in grid_params], indexing='ij')]))
pd.DataFrame(grid).to_csv(
    event + '/parameter_grid', sep='\t', index=False, float_format='%.5g')
pd.DataFrame(grid_nums, index=[0])[grid_params].to_csv(
    event + '/grid_metadata', sep='\t', index=False)
