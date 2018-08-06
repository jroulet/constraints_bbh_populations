"""Compute templates around an event and match them to the data to get 
their log-likelihood, marginalized over phase, amplitude and time.
Plot the quality of interpolating the likelihood on random values.

* Notation: 
    Capital letters refer to frequency domain, e.g.:
    H == htilde(f), s == s(t).
    
* Usage: 
    python compute_SNR.py <run_directory/>
"""

from pycbc import waveform
from pycbc.psd.estimate import welch
from pycbc.types.timeseries import TimeSeries
import os
import sys
from multiprocessing import Pool, cpu_count
import numpy as np
from numpy import sqrt, pi, exp, log, fft
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.signal import tukey
from scipy.interpolate import interp1d
from collections import OrderedDict
import time
start_time = time.time()

try:
    event = sys.argv[1].replace('/', '')
except IndexError:
    sys.exit('Usage:\npython compute_SNR.py <run_directory/>')

approximants = [#'IMRPhenomD',
                'SEOBNRv4_ROM',
                ]

os.system('mkdir -p {}/figures'.format(event))
with open('{0}/{0}.dat'.format(event)) as f:
    detectors = f.readline().strip().split()
LIGO = ['H1', 'L1']
nonLIGO = set(detectors) - set(LIGO)

# User-defined quantities:
T = 8. if event in ['GW151226', 'GW170608'] else 4.
fs = 4096.
dt = 1 / fs
df = 1 / T
f_lower = 10 if event != 'GW170608' else 20  # Excess noise for GW170608
N = int(T * fs)
times = np.linspace(0, T, num=N, endpoint=False)  # [s]

freqs = np.linspace(0, fs/2, fs/2/df + 1)
with open('{}/freqs.dat'.format(event), 'w') as freq_file:
    freq_file.write('\n'.join(['{:.6g}'.format(x) for x in freqs]))

# Load data, compute PSD:
S, WD_conj, PSD = {}, {}, {}
for i, det in enumerate(detectors):
    s = np.loadtxt('{0}/{0}.dat'.format(event), 
                   usecols=i, skiprows=1)  # GW data
    psd = welch(TimeSeries(s, delta_t=dt), avg_method='median-mean')
    psd_freqs = psd.sample_frequencies
    PSD[det] = np.interp(freqs, psd_freqs, psd)  # [1/Hz]
    if event == 'GW170608':
    	PSD['H1'][freqs < 30] = np.inf
    # Crop T seconds of data near the event:
    s = s[len(s)//2 - 3*N//4 : len(s)//2 + N//4]
    S[det] = fft.rfft(s * tukey(len(s), alpha=1./8)) / fs  # [1/Hz] 
    WD_conj[det] = S[det].conjugate() / PSD[det]  # Withened data

with open(event + '/grid_metadata') as f:
    grid_params = f.readline().strip().split()
# Some consistency checks:
with open(event + '/parameter_grid') as f:
    param_key = f.readline().strip().split()
assert not any('logL' in k for k in param_key), \
    'Likelihood already computed; run parameter_config.py again to clear.'
assert grid_params == param_key, \
    'The parameter_grid key does not match grid_metadata.'

# Load the tabulated function I(|z|) [arxiv 1806.10610 Eq. (12)]
zp, logIp = np.loadtxt('logI.dat', skiprows=1, unpack=True)
logI = interp1d(zp, logIp, kind='cubic', assume_sorted=True)

upsample_factor = 8 if event != 'GW170608' else 16

#  Prescriptions for assigning (chi1, chi2) given chi_eff:
def equal_chi1_chi2_chieff(q, chi_eff):
    return chi_eff, chi_eff
def maximal_chi1(q, chi_eff):
    chi1 = np.maximum(1, (1+q)*chi_eff + q)
    return chi1, (1+q)*chi_eff + q
def maximal_chi2(q, chi_eff):
    chi2 = np.maximum(1, ((1+q)*chi_eff + 1) / q)
    return (1+q)*chi_eff - q*chi2, chi2
spins_prescription = equal_chi1_chi2_chieff  # Pick one from the above

def m1_of_Mchirp_q(Mchirp, q):
    return Mchirp * (1 + q)**.2 / q**.6
def m2_of_Mchirp_q(Mchirp, q):
    return Mchirp * (1 + 1/q)**.2 * q**.6
def delay_prior(delay):  # [arxiv 1707.06101]
    return (1 - (delay/10.65e-3)**2) * (np.abs(delay)<10.012e-3)
def compute_likelihood(arg):
    """ Returns the multi-detector log-likelihood marginalized over
    (t0, amplitude, phase) for a given set of physical parameters.

    Compute a waveform template. For each detector, get the complex
    matched filter output z(t0), upsample. 
    For LIGO, two likelihoods are computed by combining zH and zL 
    coherently (assuming co-orientation) or incoherently (discarding 
    orientation information). Additional detectors (Virgo) are always 
    combined incoherently.
    Get the likelihood L(t0) marginalized analitically over phase and1
    amplitude interpolating the analytical (quadrature) solution.
    Marginalize numerically over t0, and for LIGO over t_delay too.
    Returns a dictionary with the log-likelihood computed coherently and
    incoherently.
    """
    approximant, p_dict = arg
    m1 = m1_of_Mchirp_q(p_dict['M_chirp'], p_dict['q'])
    m2 = m2_of_Mchirp_q(p_dict['M_chirp'], p_dict['q'])
    chi1, chi2 = spins_prescription(p_dict['q'], p_dict['chi_eff'])
    H = np.array(waveform.get_fd_waveform(
        mass1=m1,
        mass2=m2,
        spin1z=chi1,
        spin2z=chi2,
        approximant=approximant, f_lower=f_lower, f_final=fs/2, delta_f=df)[0])

    abs_z_1det_thresh = 6  # Ignore parameters that yield SNR below this
    
    dt_ups = dt / upsample_factor
    t_ups = np.arange(times[0], times[-1], dt_ups)
    logL, z, abs_z, abs_z_above_thresh, hh = {}, {}, {}, {}, {}
    for det in detectors:
        hh[det] = 4*(abs(H)**2 / PSD[det]).sum()*df  # < h | h >
        z_ = 4/sqrt(hh[det]) * fs*fft.ifft(WD_conj[det]*H, N)
        abs_z_ = np.abs(z_)
        if det in LIGO and abs_z_.max() < abs_z_1det_thresh:
            return {'logL_coherent': -99, 'logL_incoherent': -99}  # Abort
        # abs, arg vary more smoothly than Re, Im; better for interpolation:
        abs_z[det] = interp1d(times, abs_z_, kind='cubic',
                              assume_sorted=True)(t_ups)
        arg_z_ups = interp1d(times, np.unwrap(np.angle(z_)),
                             kind='cubic', assume_sorted=True)(t_ups)
        z[det] = abs_z[det] * exp(1j*arg_z_ups)
    
    # Compute the LIGO likelihood coherently, using that H1 and L1 are aligned:
    for det in LIGO:
        abs_z_above_thresh[det] = abs_z[det] > abs_z_1det_thresh
    hh['LIGO'] = hh['H1'] + hh['L1']
    f = {det: sqrt(hh[det] / hh['LIGO']) for det in LIGO}
    i_det = {det: abs_z_above_thresh[det].nonzero()[0] for det in LIGO}
    iH, iL = np.meshgrid(i_det['H1'], i_det['L1'], indexing='ij')
    delay = dt_ups * (iL - iH)
    abs_z_LIGO = np.abs(
          f['H1'] * z['H1'][abs_z_above_thresh['H1']][:, np.newaxis]
        - f['L1'] * z['L1'][abs_z_above_thresh['L1']][np.newaxis, :])
    logL_LIGO = log((exp(logI(abs_z_LIGO))*delay_prior(delay)).sum()*dt_ups**2)

    # Incoherent in LIGO, to compare:
    for det in LIGO:
        logL[det] = log(exp(logI(abs_z[det][abs_z_above_thresh[det]])).sum()
                        * dt_ups)

    # Compute the likelihood in VIRGO if it's there, combine incoherently
    di = int(30e-3 / dt_ups)  # LIGO-VIRGO distance in units of dt_ups
    for det in nonLIGO:
        look = np.full_like(t_ups, False, dtype=bool)
        for i in i_det['H1']:
            look[i-di : i+di] = True
        logL[det] = log(exp(logI(np.maximum(0, abs_z[det][look]))).sum()*dt_ups)

    logL_coherent = logL_LIGO + sum([logL[det] for det in nonLIGO])
    logL_incoherent = sum([logL[det] for det in detectors])
    return {'logL_coherent': logL_coherent,
            'logL_incoherent': logL_incoherent,
           }

p = Pool(cpu_count())

# COMPUTE LIKELIHOOD ON THE PARAMETER GRID
# ----------------------------------------
grid = OrderedDict(zip(grid_params, np.loadtxt(
    event + '/parameter_grid', skiprows=1, unpack=True)))
n_gridpoints = len(grid.values()[0])
for approximant in approximants:
    results = p.map(
        compute_likelihood, 
        [[approximant, {par: grid[par][n] for par in grid_params}]
         for n in range(n_gridpoints)])  # Compute logL in parallel
    for coherence in ['coherent', 'incoherent']:
        grid['logL', coherence, approximant] = np.array(
            [x['logL_' + coherence] for x in results])
# Save results:
joined_keys = [key if isinstance(key, basestring) else '_'.join(key) 
               for key in grid]
# For log scale, absolute rather than relative precision is relevant:
formats = ['{:.3f}' if 'logL' in key else '{:.5g}' for key in joined_keys]
with open(event + '/parameter_grid', 'w') as outfile:
    outfile.write('\t'.join(joined_keys) + '\n')
    outfile.write('\n'.join(['\t'.join([fmt.format(x) 
                                        for fmt, x in zip(formats, row)]) 
                             for row in np.transpose(grid.values())]))
print 'Time to compute logL: {:.1f} minutes'.format((time.time()-start_time)/60)

# TEST QUALITY OF INTERPOLATION WITH RANDOM PARAMETERS
# ----------------------------------------------------
# Load LIGO reported parameters
with open('{0}/{0}_LIGO_parameters'.format(event)) as f:
    lp_col_key =  f.readline().strip().split()  # lp stands for ligo_params
    lp_row_key = [l.strip().split()[0] for l in f.readlines()]
lp_cols = np.loadtxt('{0}/{0}_LIGO_parameters'.format(event), unpack=True,
                     skiprows=1, usecols=range(1, len(lp_col_key) + 1))
ligo_params = {par: dict(zip(lp_row_key, col)) 
               for par, col in zip(lp_col_key, lp_cols)}
min_val = {par: grid[par].min() for par in grid_params}
max_val = {par: grid[par].max() for par in grid_params}

# Generate random points in the LIGO credible bounds:
n_random_points = 512
mean = [ligo_params[par]['Overall'] for par in grid_params]
cov = np.diag([((  ligo_params[par]['Overall_errp'] 
                 + ligo_params[par]['Overall_errm']) / 2) ** 2 
               for par in grid_params])
random_points = np.random.multivariate_normal(mean, cov, 4*n_random_points)
random_points = [point for point in random_points 
                 if all(val > min_val[par] and val < max_val[par] 
                        for val, par in zip(point, grid_params))
                ][:n_random_points]
assert all(val > min_val[par] and val < max_val[par] 
           for point in random_points for val, par in zip(point, grid_params))
n_random_points = len(random_points)  # Just in case a lot were outliers
random_points = dict(zip(grid_params, np.transpose(random_points)))

# Compute the likelihood on the random points, compare to interpolated:
for approximant in approximants:
    results = p.map(
       compute_likelihood, 
       [[approximant, {par: random_points[par][n] for par in grid_params}]
        for n in range(n_random_points)])
    for coherence in ['coherent', 'incoherent']:
        random_points['logL_true', coherence, approximant] = np.array(
            [x['logL_' + coherence] for x in results])
    # Define the 1d grids and reconstruct the full meshgrid
    with open(event + '/grid_metadata') as grid_metadata:
        grid_params = grid_metadata.readline().strip().split()
        grid_num = OrderedDict(zip(
            grid_params, 
            [int(x) for x in grid_metadata.readline().strip().split()]))
    grid_nums = np.array(list(grid_num.values()))
    grid_1d, meshgrid = {}, {}
    for par in grid_params:
        grid_1d[par] = np.linspace(min_val[par], max_val[par], grid_num[par])
        meshgrid[par] = np.reshape(grid[par], (grid_nums))
    assert np.allclose(
        [meshgrid[par] for par in grid_params], 
        np.meshgrid(*[grid_1d[par] for par in grid_params], indexing='ij'),
        atol=1e-2), 'There was a problem reconstructing the parameter_grid'
    # Interpolate logL on the random points using the grid data
    for coherence in ['coherent', 'incoherent']:
        random_points['logL_interp', coherence, approximant] \
            = ndimage.map_coordinates(
                np.reshape(grid['logL', coherence, approximant], (grid_nums)),
                [(grid_num[par]-1) / (max_val[par]-min_val[par])
                 * (random_points[par]-min_val[par]) for par in grid_params])
        max_error = max([abs(true - interp) for true, interp in zip(
            random_points['logL_true', coherence, approximant],
            random_points['logL_interp', coherence, approximant])])
        mean_error = sqrt(np.mean([(true-interp) ** 2 for true, interp in zip(
            random_points['logL_true', coherence, approximant], 
            random_points['logL_interp', coherence, approximant])]))
        print event, coherence, approximant
        print 'Relative error in likelihood interpolation:'
        print ' max:  {:.3f}\n mean: {:.3f}'.format(max_error, mean_error)

# Plot true vs interpolated L:
fig = plt.figure()
fig.add_subplot(111, aspect='equal')
for approximant in approximants:
    for coherence in ['coherent', 'incoherent']:
        plt.plot(exp(random_points['logL_true', coherence, approximant]), 
                 exp(random_points['logL_interp', coherence, approximant]),
                 '.', ls='None', label='{}, {}'.format(coherence, approximant))
lim = plt.gca().get_xlim()
plt.plot(lim, lim, c='lightgrey', zorder=0)  # y=x line
plt.xlim(lim)
plt.ylim(lim)
plt.xlabel(r'True $\mathcal{L}$')
plt.ylabel(r'Interpolated $\mathcal{L}$')
plt.legend()
plt.savefig(event + '/figures/L_true_vs_interp.pdf', bbox_inches='tight')

print 'Total time: {:.1f} minutes'.format((time.time() - start_time) / 60)
