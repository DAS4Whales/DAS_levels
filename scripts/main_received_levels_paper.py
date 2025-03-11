import das4whales as dw
import scipy.signal as sp
import numpy as np
from utils.data_handle import find_next_file, load_data, select_channels_m_to_nb
from utils.dsp import (time_compensation, calc_g_gauge_length_db, calc_g_coupling_db,
                       calc_transmission_loss, calc_parameters, calc_received_strain_level_db, process_noise)
from utils.plot import plot_save_TOA_compensated, get_ticks_dist
import pandas as pd
import os
import ast  # For safely parsing list and tuple strings
import matplotlib.pyplot as plt

# Select the dataset to work with 'OOISouthC1' or 'Svalbard' or 'MedSea'
dataset = 'MedSea'

# Plot style
fontsize = 13
fontsize_vals = 13
fontsize_legend =13

# Analysis parameters
analysis_parameters = {
    'fs_analysis': 200, # Hz,
    'window_toa_compensation_s': 10, # (s) Window for the TOA compensation
    'window_response_s': 2, # (s) Widow for response analysis
    'sound_speed': 1490, # m/s
    'whale_depth_m': 20, # Whale vocalizing depth (m)
    'whale_freq': 20, # Hz
    'whale_SL_dB': 189, # dB re 1uPa @1m
    'alpha_cable_fiber': 0.8, #
    'bulk_water': 2.29 * 10 ** 9,
}

# Define path to folder with data
data_dir = '../data'

# Read csv file dataset_config.csv to get the relevant information dor each dataset
all_dataset_configs = pd.read_csv(os.path.join(data_dir, 'dataset_config.csv'))
print(f'Columns of the dataset configuration csv: {list(all_dataset_configs.columns)}')

# Read column selected_channels_m and sensitivity correctly
all_dataset_configs['selected_channels_m'] = all_dataset_configs['selected_channels_m'].apply(
    lambda x: ast.literal_eval(x) if pd.notna(x) else None
)
all_dataset_configs['sensitivity'] = all_dataset_configs['sensitivity'].apply(
    lambda x: ast.literal_eval(x) if pd.notna(x) else None
)

# Ensure the 'Dataset' column is set as the index
all_dataset_configs.set_index('dataset', inplace=True)
dataset_config = all_dataset_configs.loc[dataset]
analysis_parameters['H_multiple'] = dataset_config['H_multiple']

# Load annotation file to get whale_apex_m, whale_offset_m and start_time_s from it
annotations = dw.data_handle.load_annotation_csv(os.path.join(data_dir,dataset_config['annotation_file']))

# Store the information on the whale's position in a dictionary
whale_position = {
    'apex': annotations['apex'][0],
    'offset': annotations['offset'][0],
    'side': annotations['whale_side'][0],
    'depth': analysis_parameters['whale_depth_m']
     }

# Store file and start time
file = annotations['file_name'][0]
start_time_s = annotations['start_time'][0]

# Store the information on the noise file
noise = {
    'file': dataset_config['das_noise_file'],
    'start_time_s': dataset_config['noise_start_time_s']
}

# Part 1: Time compensate the labeled data
# 1) Load and condition the das data
print(f'Processing ... {os.path.basename(file)}')

if dataset_config['file_duration_s']-start_time_s < analysis_parameters['window_toa_compensation_s']:
    file_list = find_next_file(os.path.join(data_dir,file))
else:
    file_list = os.path.join(data_dir,file)

# Load the labeled file and the next to avoid edge effects in the time compensation.
tr, time, dist, fileBeginTimeUTC, metadata, selected_channels_m = load_data(dataset_config['interrogator'],
                                                                            file_list,
                                                                            dataset_config['selected_channels_m'],
                                                                            sensitivity=dataset_config['sensitivity'])

fs_original, dx, nx, ns, gauge_length = metadata["fs"], metadata["dx"], metadata["nx"], metadata["ns"], \
    metadata["GL"]
selected_channels = select_channels_m_to_nb(selected_channels_m, metadata)

# Downsample
tr, fs, time = dw.dsp.resample(tr, fs_original, analysis_parameters['fs_analysis'])

# Filter
# Create the f-k filter
fk_filter = dw.dsp.hybrid_ninf_filter_design((tr.shape[0], tr.shape[1]), selected_channels, dx, fs,
                                             cs_min=1300, cp_min=1450, cp_max=3300, cs_max=3450, fmin=14,
                                             fmax=30,
                                             display_filter=False)

# Apply the f-k filter to the data, returns spatio-temporal strain matrix
trf_fk = dw.dsp.fk_filter_sparsefilt(tr, fk_filter, tapering=False)

# 2) TOA compensation -----
# Get DAS position
position = dw.data_handle.get_cable_lat_lon_depth(
    os.path.join(data_dir,dataset_config['position_file']), selected_channels)

# Theoretical TOA
toa_th = dw.loc.calc_theory_toa(position, whale_position, dist)

# TOA compensation to get the response
trf_fk_aligned = time_compensation(trf_fk, analysis_parameters, toa_th, start_time_s)

# Save the time-compensated image
filename = 'TOA_compensated_' + fileBeginTimeUTC.strftime("%Y%m%d-%H%M%S") + '.png'
plot_save_TOA_compensated(sp.hilbert(trf_fk_aligned, axis=1), dist,
                          analysis_parameters['window_toa_compensation_s'],
                          file_begin_time_utc=fileBeginTimeUTC,
                          fig_size=(12, 10), v_min=0, v_max=0.4, ticks_distance_km=get_ticks_dist(selected_channels_m),
                          Save=[True, filename])
print(f"figure {filename} saved")

#  Response on the fiber
trf_fk_align_call = trf_fk_aligned[:,
                    int(analysis_parameters['window_toa_compensation_s'] * fs):
                    int((analysis_parameters['window_toa_compensation_s']
                         + analysis_parameters[ 'window_response_s']) * fs)]

# 3) Process the noise the same way as the labeled signal
trf_fk_align_noise = process_noise(noise, data_dir, toa_th, dataset_config, analysis_parameters,
                                   time_dim = np.shape(trf_fk)[1], plot=True)

# Part 2: Response -----------
LW = gauge_length

# 1) Calculate the distances and angles
theta, cos_phi, r, ind_whale_apex = calc_parameters(position, analysis_parameters, whale_position, dist)

# 2)  Received level at the cable
# transmisison loss
transmisison_loss_db_dict = calc_transmission_loss(position, cos_phi, r, analysis_parameters, ind_whale_apex)

# Conversion strain to Pa
g_strain_to_pa_db = 20 * np.log10(analysis_parameters['bulk_water'])

# Coupling between cable and fiber
g_coupling_db = calc_g_coupling_db(theta, analysis_parameters)

# Effect of the gauge length
g_gauge_length_db = calc_g_gauge_length_db(theta, analysis_parameters, gauge_length, LW)

# 3) Receive level on DAS (strain)
received_level_db_dict = calc_received_strain_level_db(analysis_parameters, transmisison_loss_db_dict, g_gauge_length_db,
                                                  g_coupling_db, g_strain_to_pa_db)
# DAS Response
DAS_response = g_gauge_length_db + g_coupling_db#- attenuation_DAS * dist / 1000

# 4) Receive level on DAS (dB Pa)
safe_val = 1e-10
RL_DAS_Pa_dB = (20 * np.log10(np.clip(np.mean(abs(trf_fk_align_call), axis=1) / 1e-6, safe_val, None))
                - DAS_response + g_strain_to_pa_db)



# Plot step by step  -------------------------------------
# Create a figure with 3 subplots
if dataset == 'OOISouthC1':
    width = 6.1
else:
    width = 5.45

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(width, 10), facecolor='none')

# Subplot 1: Depth vs dist
ax1.plot(dist/ 1000, abs(np.array(position['depth'])), color='orangered')
ax1.set_ylabel('Depth (m)', color='k', fontsize = fontsize)
ax1.invert_yaxis()
ax1.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

# Move tick labels to right side, inside plot
ax1.yaxis.tick_right()  # Move ticks to right side
ax1.tick_params(axis='y', direction='in', pad=-5, labelsize=fontsize_vals)  # Move labels inside
ax1.yaxis.set_label_position('left')  # Keep ylabel on left

# Right-align the tick labels
for tick in ax1.yaxis.get_major_ticks():
    tick.label2.set_horizontalalignment('right')  # Right-align the labels

ax1.grid()

# Subplot 2: Acoustic received levels on the fiber
# Transmission loss and Lloyds Mirror effect
ax2.plot(dist / 1000,
         analysis_parameters['whale_SL_dB'] + transmisison_loss_db_dict['spherical'],
         'steelblue', label=r'$RL_{dB \; re. 1\mu Pa} \; (r_T \to \infty)$')
ax2.plot(dist / 1000,
         analysis_parameters['whale_SL_dB'] + transmisison_loss_db_dict['cylindrical'] - 30,
         'steelblue', linestyle='dashed',
         label=r'$RL_{dB \; re. 1\mu Pa} \; (r_T \to 1) - 30 \;$ dB')
ax2.plot(dist / 1000,
         analysis_parameters['whale_SL_dB'] + transmisison_loss_db_dict['regime_shift'],
         'k', label=r'$RL_{dB \; re. 1\mu Pa}$')

# Labels
ax2.tick_params( labelsize=fontsize_vals)
ax2.set_ylabel('Rec. level (dB re. 1$\mu$Pa)', fontsize = fontsize)
ax2.yaxis.set_label_coords(-0.11, 0.5)

ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax2.set_xlabel('')  # Remove x-axis label for ax2

# Limits
ax2.set_ylim([95, 153])  # For the global plot
ax2.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax2.grid()

# Subplot 3: DAS Response (strain)
# Coupling
ax3.plot(dist / 1000, g_coupling_db, 'maroon', label='$G_{\kappa}$')
# Gauge length
ax3.plot(dist / 1000, g_gauge_length_db, 'indianred', linestyle='dashdot', label='$G_{GL}$')
# Combined
ax3.plot(dist / 1000, DAS_response, 'k', label=r'$H_{DAS_{\varepsilon}}$')

# Labels
ax3.set_ylabel(r'DAS resp. (dB re. 1$\mu\varepsilon$)', fontsize = fontsize)
ax3.yaxis.set_label_coords(-0.11, 0.5)
ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax3.set_xlabel('')  # Remove x-axis label for ax3
ax3.tick_params(labelsize=fontsize_vals)

# Limits
ax3.set_ylim([-25, 5])  # For the global plot
ax3.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax3.grid()

# Subplot 4: Total RL reference in strain
# Plot the time-compensated noise
ax4.plot(dist / 1000, 20 * np.log10(np.clip(np.mean(abs(trf_fk_align_noise), axis=1) / 1e-6, safe_val, None)) ,
                     label=r'$NL_{dB \; re. 1\mu \varepsilon} \;$ Measured', color='slategrey', alpha = 0.8, linewidth = 0.5)
# Plot the time-compensated call
ax4.plot(dist / 1000, 20 * np.log10(np.clip(np.mean(abs(trf_fk_align_call), axis=1) / 1e-6, safe_val, None)) ,
         label=r'$RL_{dB \; re. 1\mu \varepsilon} \;$ Measured', color='orangered')
# Plot the simulations with regime shift
ax4.plot(dist / 1000, received_level_db_dict['regime_shift'], 'k',
         label=r'$RL_{dB \; re. 1\mu \varepsilon} \;$ Simulated, $S_{Pa \to \varepsilon} = 187.2$ dB')


# Labels
ax4.tick_params(labelsize=fontsize_vals)
ax4.set_ylabel(r'Rec. level (dB re. 1$\mu\varepsilon$)', fontsize = fontsize)
ax4.yaxis.set_label_coords(-0.11, 0.5)
ax4.set_xlabel('Distance (km)', fontsize = fontsize)

ax4.set_ylim([-105, -50])  # For the global plot
ax4.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax4.grid()

if dataset  == 'Svalbard':
    ax1.set_xlim([50 , dist[-1] / 1000])
    ax2.set_xlim([50, dist[-1] / 1000])
    ax3.set_xlim([50, dist[-1] / 1000])
    ax4.set_xlim([50, dist[-1] / 1000])
    ax2.legend(loc='upper right', fontsize = fontsize_legend)
    ax3.legend(loc='lower right', fontsize = fontsize_legend)
    ax4.legend(loc='upper right', fontsize = fontsize_legend)
    ax1.set_ylim([350, 0])
    manual_ticks = [50, 100, 150, 200, 250, 300]  # Adjust these values as needed
elif dataset == 'OOISouthC1':
    ax1.set_ylim([650, 0])
    manual_ticks = [100, 200, 300, 400, 500, 600]  # Adjust these values as needed
elif dataset == 'MedSea':
    ax1.set_ylim([2500, 0])
    manual_ticks = [500, 1000, 1500, 2000]  # Adjust these values as needed

ax1.set_yticks(manual_ticks)

if dataset != 'OOISouthC1':
    ax1.set_ylabel('')  # Remove x-axis label for ax2
    ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax2.set_ylabel('')  # Remove x-axis label for ax2
    ax3.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax3.set_ylabel('')  # Remove x-axis label for ax2
    ax4.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax4.set_ylabel('')  # Remove x-axis label for ax2

# Adjust layout and save
plt.tight_layout()
filename = os.path.join(f'{dataset}_Receive_Level_combined-{os.path.basename(file)}'+ '.png')
plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=180)
plt.show(block=False)
plt.clf()
plt.close(fig)
print(f"figure {filename} saved")