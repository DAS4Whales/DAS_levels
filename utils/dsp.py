import das4whales as dw
import scipy.signal as sp
import numpy as np
from utils.plot import plot_save_TOA_compensated, get_ticks_dist
import os

def value_mod_1(x):
    """
    Normalizes the input value to be within the range [-1, 1].

    Parameters:
    - x: Input value (can be scalar or array).

    Returns:
    - Normalized value(s) wrapped to the range [-1, 1].
    """
    return ((x + 1) % 2) - 1

def time_compensation(trf_fk, analysis_parameters, toa_th, start_time_s, plot=False):

    # fs
    fs = analysis_parameters['fs_analysis']

    # Check if the start time duration + analysis_parameters['window_toa_compensation_s'] does not exceed file duration
    # Get the signal portion that correspond to the analysis window
    trf_fk_call = trf_fk[:,
                  int(np.floor(start_time_s * fs)):int(
                      np.floor((start_time_s + analysis_parameters['window_toa_compensation_s']) * fs))]

    # Now we want to create a matrix where we will do the compensation:
    trf_fk_aligned = np.zeros([int(np.shape(trf_fk)[0]), int(2 * analysis_parameters['window_toa_compensation_s'] * fs)])

    for dd in range(len(toa_th)):
        # Calculate the number of zeros to prepend and append
        if int(np.round(analysis_parameters['window_toa_compensation_s'] * fs - toa_th[dd] * fs)) > 0:
            zeros_prepend = np.zeros([int(np.round(analysis_parameters['window_toa_compensation_s'] * fs - toa_th[dd] * fs))])
        else:
            zeros_prepend = np.zeros([0])
        zeros_append = np.zeros([int(toa_th[dd] * fs)])

        # Ensure we have the correct size before concatenating
        if int(2 * analysis_parameters['window_toa_compensation_s'] * fs) - (
                len(zeros_prepend) + len(trf_fk_call[dd, :]) + len(zeros_append)) != 0:
            zeros_append = np.zeros([int(toa_th[dd] * fs) + (int(2 * analysis_parameters['window_toa_compensation_s'] * fs) - (
                    len(zeros_prepend) + len(trf_fk_call[dd, :]) + len(zeros_append)))])

        # Compensate for the selected row
        trf_fk_aligned[dd, :] = np.concatenate((zeros_prepend, trf_fk_call[dd, :], zeros_append), axis=0)

    return trf_fk_aligned