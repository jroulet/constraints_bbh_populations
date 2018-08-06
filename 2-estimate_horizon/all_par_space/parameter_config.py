"""Write a configuration file `parameter_grid` with parameter values on
which `compute_SNR.py` will generate templates.
"""

import numpy as np
import pandas as pd
import sys
from collections import OrderedDict

grid_nums = OrderedDict([('M_chirp', 128),
                         ('q'      , 16),
                         ('chi_eff', 64)
                        ])  # Number of points along each dimension
bounds = {'M_chirp': [4.3, 100],
          'q':       [.03, 1],
          'chi_eff': [-1, 1]
         }
grid_params = list(grid_nums.keys())
grid_1d = OrderedDict()
for par in grid_params:
    grid_1d[par] = np.linspace(*bounds[par], grid_nums[par])

grid = OrderedDict(zip(grid_params, [x.flatten() for x in np.meshgrid(
    *[grid_1d[par] for par in grid_params], indexing='ij')]))
# Write grid:
pd.DataFrame(grid).to_csv(
    'parameter_grid', sep='\t', index=False, float_format='%.5g')
pd.DataFrame(grid_nums, index=[0])[grid_params].to_csv(
    'grid_metadata', sep='\t', index=False)