# Supplementary material to "Constraints on Binary Black Hole Populations from LIGO-Virgo Detections"

This repository contains the code used to reanalyze the LIGO-Virgo data and do the population inference described in arxiv.org/abs/1806.10610

It is free to use.

## Purpose

1. Estimate parameter likelihood of BBH merger events using a grid in parameter space (chirp mass, mass ratio, effective spin).

2. Estimate horizon distance dependence on parameters to account for selection effects.

3. Infer population parameters describing the BBH merger rate dependence on the BBH parameters.

## Requisites

pycbc, jupyter

## Usage
#### 1-estimate_parameters
* Write a file `<event>.dat` whose columns are the detector strains sampled at 4096 Hz, and header `H1   L1  ...`
* Write a file `<event>_ligo_parameters` with the parameters LIGO reported for that event.  
*Tip:* Edit and run `source2detector_frame.py` if LIGO only reported source-frame quantities.
* Run `python config_parameters.py <event>` for the desired event (generates grid in parameter space). Alternatively, run `bash all_config` to do all events at once.
* Run `python compute_likelihood.py <event>` within the pycbc environment (generates templates on the grid and matches them to the data). Alternatively, run `bash all_compute` to do all events at once.
* Optional: analyze the event data with `estimate_parameters.ipynb`

#### 2-estimate_horizon
* Run `config_parameters.py` (generates grid in parameter space).
* Run `compute_SNR_1Mpc.py all_par_space/` with pycbc (generates templates on the grid and stores the expected SNR of an optimally-aligned source at 1 Mpc, given a reference amplitude spectral density (ASD)).

#### 3-analyze_population
* Run `analyze_population.ipynb`


## File structure

    .
    ├── 1-estimate_parameters/				# Single-event analysis
    │   ├── all_clear.sh					# Clear all the runs (run 1st)
    │   ├── all_compute.sh					# Run compute_likelihood.py on all the events (run 3rd)
    │   ├── all_config.sh					# Run parameter_config.py on all the events (run 2nd)
    │   ├── compute_likelihood.log2
    │   ├── compute_likelihood.py			# Compute parameter likelihood over parameter_grid
    │   ├── estimate_parameters.ipynb		# Analyze the single-event likelihoods
    │   ├── logI.dat						# Tabulated function log(I(|z|))
    │   ├── logI.nb							# Generate logI.dat
    │   ├── parameter_config.py				# Generate parameter_grid and grid_metadata
    │   ├── figures/						# Multi-event figures
    │   ├── <event1>/
    │   │   ├── figures/					# Single-event figures
    │   │   ├── freqs.dat					# Template frequencies (Hz)
    │   │   ├── grid_metadata				# Dimensions of the parameter grid
    │   │   ├── <event1>.dat				# Strain data
    │   │   ├── <event1>_LIGO_parameters	# Parameter values reported by the LVC
    │   │   └── parameter_grid				# Parameter values on which to estimate the likelihood
    │   └── .../
    │
    ├── 2-estimate_horizon/
    │   ├── ASD_average.dat					# A reference amplitude spectral density
    │   ├── ASD_zdhp.dat
    │   ├── compute_PSD_average.py			# Generate ASD_average.dat
    │   ├── compute_SNR_1Mpc.py				# Compute the SNR over parameter_grid
    │   └── all_par_space/
    │       ├── figures/
    │       ├── freqs.dat
    │       ├── grid_metadata
    │       ├── parameter_config.py
    │       └── parameter_grid
    │
    ├── 3-analyze_population
    │   ├── analyze_population.ipynb		# Do population inference
    │   ├── P_greater_than_w.dat			# Tabulated cumulative of angular factors, P(w > w*)
    │   ├── P_greater_than_w.ipynb			# Generate P_greater_than_w.dat
    │   └── figures/
    │
    └── README.md							# This file



***
## Downloading data
Download the event data in txt format from https://losc.ligo.org/events and merge the detector files, e.g. with

    `paste H1_file L1_file ... > event.dat`
Then edit the headers so they read `H1    L1    ...`

## Pycbc
Install e.g. with

    docker pull pycbc/pycbc-el7:latest
Run e.g. with

    docker run -e DISPLAY=${DISPLAY} -v /tmp/.X11-unix:/tmp/.X11-unix/ -v ${HOME}/<path/to>/constraints_bbh_populations:/home/pycbc/constraints_bbh_populations -v ${HOME}/.ssh:/home/pycbc/.ssh -it pycbc/pycbc-el7:latest /bin/bash -l
