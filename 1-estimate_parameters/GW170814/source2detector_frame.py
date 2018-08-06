import pandas as pd
from collections import OrderedDict

event = 'GW170814'
# _s : source frame
m1_s          = 30.5
m1_errm_s     = 3
m1_errp_s     = 5.7
m2_s          = 25.3
m2_errm_s     = 4.2
m2_errp_s     = 2.8
Mchirp_s      = 24.1
Mchirp_errm_s = 1.1
Mchirp_errp_s = 1.4
z             = .11
z_errm        = .04
z_errp        = .03
chi_eff       = .06
chi_eff_errp  = .12
chi_eff_errm  = .12 

Mchirp = (1 + z) * Mchirp_s
Mchirp_errm = (1 + z) * (Mchirp_errm_s**2 - (Mchirp_s*z_errp))**.5
Mchirp_errp = (1 + z) * (Mchirp_errp_s**2 - (Mchirp_s*z_errm))**.5

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