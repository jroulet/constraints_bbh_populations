"""Compute a nominal PSD to use for sensitive-volume estimation, by averaging
all the events' PSDs as described in arxiv 1806.10610.
Write the ASD=sqrt(PSD) to file.
"""
from pycbc import waveform
from pycbc.psd.estimate import welch
from pycbc.types.timeseries import TimeSeries
import sys
import numpy as np
from numpy import sqrt, pi, exp, log, fft
import matplotlib.pyplot as plt
from scipy.signal import tukey
from scipy.interpolate import interp1d

events = ['GW150914',
          'GW151226',
          'LVT151012',
          'GW170104',
          'GW170608',
          'GW170814'
         ]

LIGO = ['H1', 'L1']

# User-defined quantities:
T = 8.
fs = 4096.
dt = 1 / fs
df = 1 / T
f_lower = 10
N = int(T * fs)
times = np.linspace(0, T, num=N, endpoint=False)  # [s]

freqs = np.linspace(0, fs/2, fs/2/df + 1)

# Load data, compute PSD:
PSD = {}
for event in events:
    for i, det in enumerate(LIGO):
        s = np.loadtxt('../1-estimate_parameters/{0}/{0}.dat'.format(event), 
                       usecols=i, skiprows=1)  # GW data
        psd = welch(TimeSeries(s, delta_t=dt), avg_method='median-mean')
        psd_freqs = psd.sample_frequencies
        PSD[event, det] = np.interp(freqs, psd_freqs, psd)  # [1/Hz]

    PSD[event, 'total'] = (PSD[event, 'H1']**-2 + PSD[event, 'L1']**-2)**-.5
PSD['average'] = 1/np.mean([1/PSD[event, 'total'] for event in events], axis=0)

ASD = {x: sqrt(PSD[x]) for x in PSD}
with open('ASD_average.dat', 'w') as outfile:
    outfile.write('\n'.join('{}\t{}'.format(f, a) 
                            for f, a in zip(freqs, ASD['average'])))

fmin = 20
fmax = 1500
imin = next(i for i, f in enumerate(freqs) if f > fmin)
imax = next(i for i, f in enumerate(freqs) if f > fmax)
for x in ASD:
    ASD[x] = ASD[x][imin:imax]
freqs = freqs[imin:imax]

plt.loglog(freqs, ASD['average'], 'k', label='Average', lw=2, zorder=10)
for event in events:
    plt.plot(freqs, ASD[event, 'total'], label=event)
plt.xlim(fmin)
plt.xlabel('Frequency (Hz)')
plt.ylabel('ASD (Hz$^{-1/2}$)')
plt.legend()
plt.savefig('ASD_average.pdf', bbox_inches='tight')
plt.show()