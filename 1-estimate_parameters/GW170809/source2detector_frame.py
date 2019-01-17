"""Reverse-engineer the detector frame values and uncertainties from the
source frame quantities LIGO reported.
"""

import pandas as pd
from collections import OrderedDict

event = 'GW170809'
# _s : source frame
m1_s          = 35.2
m1_errm_s     = 6.
m1_errp_s     = 8.3
m2_s          = 23.8
m2_errm_s     = 5.1
m2_errp_s     = 5.2
Mchirp_s      = 25.
Mchirp_errm_s = 1.6
Mchirp_errp_s = 2.1
z             = .2
z_errm        = .07
z_errp        = .05
chi_eff       = .07
chi_eff_errm  = .16
chi_eff_errp  = .16

Mchirp = (1 + z) * Mchirp_s
#Mchirp_errm = (1 + z) * (Mchirp_errm_s**2 - (Mchirp_s*z_errp))**.5
#Mchirp_errp = (1 + z) * (Mchirp_errp_s**2 - (Mchirp_s*z_errm))**.5
Mchirp_errm = 0
Mchirp_errp = 0

q = m2_s / m1_s
q_errm = ((m2_errm_s/m1_s)**2 + (m2_s/m1_s**2 * m1_errp_s)**2)**.5
q_errp = ((m2_errp_s/m1_s)**2 + (m2_s/m1_s**2 * m1_errm_s)**2)**.5

p = OrderedDict((('M_chirp', [Mchirp, Mchirp_errm, Mchirp_errp]),
                ('q', [q, q_errm, q_errp]),
                ('chi_eff', [chi_eff, chi_eff_errm, chi_eff_errp])
               ))


ligo_parameters = pd.DataFrame(p,
	index=['Overall', 'Overall_errm', 'Overall_errp'])
ligo_parameters.to_csv(event + '_LIGO_parameters', 
	                   sep='\t', float_format='%.4g')