"""Reverse-engineer the detector frame values and uncertainties from the
source frame quantities LIGO reported.
"""

import pandas as pd
from collections import OrderedDict

event = 'GW170823'
# _s : source frame
m1_s          = 39.6
m1_errm_s     = 6.6
m1_errp_s     = 10
m2_s          = 29.4
m2_errm_s     = 7.1
m2_errp_s     = 6.3
Mchirp_s      = 29.3
Mchirp_errm_s = 3.2
Mchirp_errp_s = 4.2
z             = .34
z_errm        = .14
z_errp        = .13
chi_eff       = .08
chi_eff_errp  = .2
chi_eff_errm  = .22

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