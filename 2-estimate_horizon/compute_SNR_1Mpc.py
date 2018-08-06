"""Compute templates on a grid of parameters. Match filter
with themselves given some reference PSD to estimate the
SNR that the event would have at 1Mpc and optimal orientation.
Plot the quality of interpolating the SNR on random values.
*Notation:* capital letters refer to frequency domain, e.g.:
    H == htilde(f), s == s(t).
"""

from pycbc import waveform
import os
import sys
from multiprocessing import Pool, cpu_count
import numpy as np
from numpy import sqrt, pi, exp, log, fft
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.signal import tukey
from scipy.interpolate import interp1d
import scipy.ndimage as ndimage
from collections import OrderedDict
import cmath
import time
start_time = time.time()

try:
    path = sys.argv[1].replace('/', '')
except IndexError:
    sys.exit('Usage:\npython compute_SNR.py <run_directory/>')

approximants = ['IMRPhenomD', 'SEOBNRv4_ROM']
PSDs = ['ASD_zdhp',
        'ASD_average']

os.system('mkdir -p {}/figures'.format(path))

# User-defined quantities:
T = 8.
fs = 4096.
dt = 1 / fs
df = 1 / T
f_lower = 10
N = int(T * fs)
times = np.linspace(-T, 0, num=N, endpoint=False)  # [s]

freqs = np.linspace(0, fs/2, fs/2/df + 1)

fmin = 10
fmax = 1600  # Where the PSD is well calibrated
imin = np.argmax(freqs > fmin)
imax = np.argmax(freqs > fmax)

with open('{}/freqs.dat'.format(path), 'w') as freq_file:
    freq_file.write('\n'.join(['{:.6g}'.format(x) for x in freqs]))

def m1_of_Mchirp_q(Mchirp, q):
    return Mchirp * (1 + q)**.2 / q**.6
def m2_of_Mchirp_q(Mchirp, q):
    return Mchirp * (1 + 1/q)**.2 * q**.6

def compute_SNR_1Mpc(arg):
    '''Return the one-detector SNR an event of given parameters would 
    have, given the PSD and assuming D=1Mpc, optimal orientation and
    that the template is identical to the signal.
    '''
    approximant, p_dict = arg
    m1 = m1_of_Mchirp_q(p_dict['M_chirp'], p_dict['q'])
    m2 = m2_of_Mchirp_q(p_dict['M_chirp'], p_dict['q'])
    H = np.array(waveform.get_fd_waveform(
        mass1=m1,
        mass2=m2,
        spin1z=p_dict['chi_eff'],
        spin2z=p_dict['chi_eff'],  # Force chi1 == chi2 == chi_eff
        approximant=approximant, f_lower=f_lower, f_final=fs/2, delta_f=df)[0])
    SNR_1Mpc = {}
    for psd in PSDs:
        SNR_1Mpc[psd] = sqrt((4*abs(H[imin : imax])**2 / PSD[psd] * df).sum())
    return SNR_1Mpc

# Load PSD:
asd_freqs, asd_raw = {}, {}  # Original ASD frequencies
ASD, PSD = {}, {}  # Resampled freqs to match templates
for psd in PSDs:  # psd: name of the sensitivity curve
    asd_freqs[psd], asd_raw[psd] = np.loadtxt(psd + '.dat', unpack=True)
    ASD[psd] = np.interp(freqs, asd_freqs[psd], asd_raw[psd]) # [1/Hz]
    PSD[psd] = ASD[psd][imin : imax] ** 2

p = Pool(cpu_count())
with open(path + '/grid_metadata') as f:
    grid_params = f.readline().strip().split()
# Some consistency checks:
with open(path + '/parameter_grid') as f:
    param_key = f.readline().strip().split()
assert not any('SNR' in k for k in param_key), \
    'SNR already computed; run parameter_config.py again to clear.'
assert grid_params == param_key, \
    'The parameter_grid key does not match grid_metadata.'


# COMPUTE THE 1 Mpc SNR ON THE PARAMETER GRID
# -------------------------------------------
grid = OrderedDict(zip(grid_params, np.loadtxt(
    path + '/parameter_grid', skiprows=1, unpack=True)))
n_gridpoints = len(grid.values()[0])
for approximant in approximants:
    SNR_1Mpc = p.map(
        compute_SNR_1Mpc, 
        [[approximant, {par: grid[par][n] for par in grid_params}]
         for n in range(n_gridpoints)])  # Collect results over grid
    for psd in PSDs:
        grid['SNR_1Mpc', psd, approximant] = np.array([x[psd] for x in SNR_1Mpc])
# Save results:
joined_keys = [key if isinstance(key, basestring) else '_'.join(key) 
               for key in grid]
with open(path + '/parameter_grid', 'w') as outfile:
    outfile.write('\t'.join(joined_keys) + '\n')
    outfile.write('\n'.join(['\t'.join(['{:.5g}'.format(x) for x in row]) 
                             for row in np.transpose(grid.values())]))
print 'Time to compute 1Mpc SNR: {:.1f} minutes'.format(
    (time.time()-start_time)/60)


# TEST QUALITY OF INTERPOLATION WITH RANDOM PARAMETERS
# ----------------------------------------------------
min_val = {par: grid[par].min() for par in grid_params}
max_val = {par: grid[par].max() for par in grid_params}
# Generate random points uniformly over the grid extent:
n_random_points = 512
random_points = np.random.uniform(low=[min_val[par] for par in grid_params],
                                  high=[max_val[par] for par in grid_params],
                                  size=[n_random_points, len(grid_params)])
random_points = OrderedDict(zip(grid_params, np.transpose(random_points)))
for approximant in approximants:
    # Compute the SNR_1Mpc on the random points
    SNR_1Mpc = p.map(
        compute_SNR_1Mpc, 
        [[approximant, {par: random_points[par][n] for par in grid_params}]
         for n in range(n_random_points)])
    for psd in PSDs:
        random_points['SNR_1Mpc_true', psd, approximant] = np.array(
            [x[psd] for x in SNR_1Mpc])
    # Define the 1d grids and reconstruct the full meshgrid
    with open(path + '/grid_metadata') as grid_metadata:
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
    # Interpolate SNR_1Mpc on the random points using the grid data
    for psd in PSDs:
        random_points['SNR_1Mpc_interp', psd, approximant] \
            = ndimage.map_coordinates(
                np.reshape(grid['SNR_1Mpc', psd, approximant], grid_nums),
                [(grid_num[par]-1) / (max_val[par]-min_val[par])
                * (random_points[par]-min_val[par]) for par in grid_params])
        max_error = max([abs(true - interp)/true for true, interp in zip(
            random_points['SNR_1Mpc_true', psd, approximant], 
            random_points['SNR_1Mpc_interp', psd, approximant])])
        mean_error = np.mean(
            [(true-interp)/true for true, interp in zip(
                random_points['SNR_1Mpc_true', psd, approximant], 
                random_points['SNR_1Mpc_interp', psd, approximant])])
        print path, psd, approximant
        print 'Error in 1 Mpc SNR interpolation:'
        print ' max:  {:.3f}\n mean: {:.3f}'.format(max_error, mean_error)

# Plot true vs interpolated SNR_1Mpc:
for psd in PSDs:
    fig = plt.figure()
    fig.add_subplot(111, aspect='equal')
    for approximant in approximants:
        plt.plot(random_points['SNR_1Mpc_true', psd, approximant],
                 random_points['SNR_1Mpc_interp', psd, approximant],
                 '.', ls='None', label=approximant)
    lim = plt.gca().get_xlim()
    plt.plot(lim, lim, c='lightgrey', zorder=0)  # y=x line
    plt.xlim(lim)
    plt.ylim(lim)
    plt.xlabel('True 1 Mpc SNR')
    plt.ylabel('Interpolated 1 Mpc SNR')
    plt.title(psd)
    plt.legend()
    plt.savefig('{}/figures/SNR_1Mpc_true_vs_interp_{}.pdf'.format(path, psd),
                bbox_inches='tight')

print 'Total time: {:.1f} minutes'.format((time.time() - start_time) / 60)
