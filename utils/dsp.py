import das4whales as dw
import scipy.signal as sp
import numpy as np
from utils.plot import plot_save_TOA_compensated, get_ticks_dist
from utils.data_handle import load_data, select_channels_m_to_nb
import os
safe_val = 1e-12

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

def process_noise(noise, data_dir, toa_th, dataset_config, analysis_parameters, time_dim, plot=False):
    tr, time, dist, fileBeginTimeUTC, metadata, selected_channels_m = load_data(dataset_config['interrogator'],
                                                                                os.path.join(data_dir, noise['file']),
                                                                                dataset_config['selected_channels_m'],
                                                                                sensitivity=dataset_config[
                                                                                    'sensitivity'])
    fs_original, dx, nx, ns, gauge_length = metadata["fs"], metadata["dx"], metadata["nx"], metadata["ns"], \
        metadata["GL"]

    selected_channels = select_channels_m_to_nb(selected_channels_m, metadata)

    # Downsample
    tr, fs, time = dw.dsp.resample(tr, fs_original, analysis_parameters['fs_analysis'])

    # added these lines specifically for the MedSea data -- not all the files are the same length
    if dataset_config['interrogator'] == 'asn':
        if time_dim > np.shape(tr)[1]:
            tr = np.pad(tr, ((0, 0), (0, time_dim-np.shape(tr)[1])), mode='constant', constant_values=0)
        elif time_dim < np.shape(tr)[1]:
            tr = tr[:,:time_dim]

    # Filter
    # Create the f-k filter
    fk_filter = dw.dsp.hybrid_ninf_filter_design((tr.shape[0], tr.shape[1]), selected_channels, dx, fs,
                                                 cs_min=1300, cp_min=1450, cp_max=3300, cs_max=3450, fmin=14,
                                                 fmax=30,
                                                 display_filter=False)

    # Apply the f-k filter to the data, returns spatio-temporal strain matrix
    trf_fk = dw.dsp.fk_filter_sparsefilt(tr, fk_filter, tapering=False)

    # TOA compensation to get the response
    trf_fk_aligned = time_compensation(trf_fk, analysis_parameters, toa_th, noise['start_time_s'])

    if plot:
        # Save the time-compensated image
        filename = 'TOA_compensated_noise_' + fileBeginTimeUTC.strftime("%Y%m%d-%H%M%S") + '.png'
        plot_save_TOA_compensated(sp.hilbert(trf_fk_aligned, axis=1), dist,
                                  analysis_parameters['window_toa_compensation_s'],
                                  file_begin_time_utc=fileBeginTimeUTC,
                                  fig_size=(12, 10), v_min=0, v_max=0.4,
                                  ticks_distance_km=get_ticks_dist(selected_channels_m),
                                  Save=[True, filename])
        print(f"figure {filename} saved")

    #  Response on the fiber
    trf_fk_align_noise = trf_fk_aligned[:,
                         int(analysis_parameters['window_toa_compensation_s'] * fs):
                         int((analysis_parameters['window_toa_compensation_s']
                              + analysis_parameters['window_response_s']) * fs)]
    return trf_fk_align_noise

def calc_parameters(position, analysis_parameters, whale_position, dist):
    # Depth difference z
    z = np.abs(np.array(position['depth'])) - analysis_parameters['whale_depth_m']

    # Horizontal distance to the closest channel hy
    hy = whale_position['offset']

    # Distance between apex and channels, hx
    # Find the index of whale_apex_m in dist
    ind_whale_apex = np.where(dist >= whale_position['apex'])[0][0]

    # Correct the whale offset if 0, becomes depth-whale depth
    whale_offset_m_temp = int(
        np.floor(np.array(position['depth'])[ind_whale_apex] - analysis_parameters['whale_depth_m']))
    if whale_position['offset'] <= whale_offset_m_temp:
        whale_offset_m = whale_offset_m_temp

    # Calculate hx using the coordinates
    whale_pos_projected_DAS = {
        'lat': position['lat'][ind_whale_apex],
        'lon': position['lon'][ind_whale_apex],
    }
    hx = dw.spatial.calc_dist_lat_lon(whale_pos_projected_DAS, position)

    # Calculate r (distance from whale position to cable positions)
    # Get the bearing of the DAS cable around whale position
    step = 3
    DAS_bearing = dw.spatial.calc_das_section_bearing(
        position['lat'][ind_whale_apex - step],
        position['lon'][ind_whale_apex - step],
        position['lat'][ind_whale_apex + step],
        position['lon'][ind_whale_apex + step])

    # Get the whale position
    whale_position['lat'], whale_position['lon'] = dw.spatial.calc_source_position_lat_lon(
        position['lat'][ind_whale_apex],
        position['lon'][ind_whale_apex],
        whale_position['offset'],
        DAS_bearing,
        whale_position['side'])
    r = dw.spatial.calc_dist_lat_lon(whale_position, position)

    # Grazing angle on the fiber, theta
    theta = np.arccos(value_mod_1(hx / r))  # in rad

    # The vertical angle, Phi
    cos_phi = z / r
    return theta, cos_phi, r, ind_whale_apex

def calc_g_gauge_length_db(theta, analysis_parameters, gauge_length, LW):
    # Precompute common terms
    k = np.pi * analysis_parameters['whale_freq'] / analysis_parameters['sound_speed']  # Wave number/2
    k_gauge = k * gauge_length * np.cos(theta)
    k_LW = k * LW * np.cos(theta)

    # Sinc calculations, with handling for division by zero
    term_gauge = np.sin(k_gauge) / k_gauge #if k_gauge != 0 else 1
    term_LW = np.sin(k_LW) / k_LW #if k_LW != 0 else 1

    # Effect of the gauge length
    g_GL = term_gauge * term_LW

    return 20*np.log10(g_GL)

def calc_g_coupling_db(theta, analysis_parameters):
    g_coup_fiber_cable = (0.7 * np.cos(theta) ** 2 - 0.2 * analysis_parameters['alpha_cable_fiber'] * np.sin(theta) ** 2)
    return 20*np.log10(np.clip(g_coup_fiber_cable, safe_val, None))

def calc_transmission_loss(position, cos_phi, r, analysis_parameters, ind_whale_apex):
    g_lme = 2 * np.sin(2 * np.pi * analysis_parameters['whale_freq'] * analysis_parameters['whale_depth_m'] * cos_phi /
                       analysis_parameters['sound_speed'])
    g_spherical = 1 / r
    TL_water_spherical_dB = 20 * np.log10(g_lme) + 20 * np.log10(g_spherical)

    g_cylindrical = 1 / np.sqrt(r)
    TL_water_cylindrical_dB = 20 * np.log10(g_cylindrical)

    # Regime shift at 4 x the height of the water column at the whale location
    g_regime_shift = g_spherical
    TL_water_regime_shift_dB = 20 * np.log10(g_lme) + 20 * np.log10(g_regime_shift)

    # going in the negative dist
    try:
        r_negative = r[0:ind_whale_apex]
        ind_regime_shift_negative = ind_whale_apex - 1 - np.where(
            r_negative[::-1] >= analysis_parameters['H_multiple'] * position['depth'][ind_whale_apex])[0][0]
        TL_water_regime_shift_dB[0:ind_regime_shift_negative] = -10 * np.log10(
            r[ind_regime_shift_negative + 1]) + 20 * np.log10(
            g_lme[ind_regime_shift_negative + 1]) - 10 * np.log10(r[0:ind_regime_shift_negative])
    except:
        print('No regime shift')
    try:
        # going in the positive r
        r_positive = r[ind_whale_apex:]
        ind_regime_shift_positive = ind_whale_apex + \
                                    np.where(r_positive >= analysis_parameters['H_multiple'] * position['depth'][
                                        ind_whale_apex])[0][0]
        TL_water_regime_shift_dB[ind_regime_shift_positive:] = -10 * np.log10(
            r[ind_regime_shift_positive - 1]) + 20 * np.log10(
            g_lme[ind_regime_shift_positive - 1]) - 10 * np.log10(r[ind_regime_shift_positive:])
    except:
        print('No regime shift')

    tl_db_dict ={
        'spherical': TL_water_spherical_dB,
        'cylindrical': TL_water_cylindrical_dB,
        'regime_shift': TL_water_regime_shift_dB
    }
    return tl_db_dict

def calc_received_strain_level_db(analysis_parameters, tl_db_dict, g_gauge_length_db, g_coupling_db, g_strain_to_pa_db):
    RL_DAS_spherical_dB = analysis_parameters['whale_SL_dB'] + tl_db_dict[
        'spherical'] + g_gauge_length_db + g_coupling_db - g_strain_to_pa_db

    RL_DAS_cylindrical_dB = analysis_parameters['whale_SL_dB'] + tl_db_dict[
        'cylindrical'] + g_gauge_length_db + g_coupling_db - g_strain_to_pa_db

    RL_DAS_regime_shift_dB = analysis_parameters['whale_SL_dB'] + tl_db_dict[
        'regime_shift'] + g_gauge_length_db + g_coupling_db - g_strain_to_pa_db

    received_level_db = {
        'spherical': RL_DAS_spherical_dB,
        'cylindrical': RL_DAS_cylindrical_dB,
        'regime_shift': RL_DAS_regime_shift_dB
    }

    return received_level_db
