import das4whales as dw
import scipy.signal as sp
import numpy as np
from utils.data_handle import find_next_file, load_data
from utils.dsp import time_compensation, value_mod_1
from utils.plot import plot_save_TOA_compensated, get_ticks_dist
import pandas as pd
import os
import ast  # For safely parsing list and tuple strings
import matplotlib.pyplot as plt

fontsize = 14
fontsize_vals = 13
fontsize_legend =13

# Analysis parameters
analysis_parameters = {
    'fs_analysis': 200, # Hz,
    'window_toa_compensation_s': 10, # (s) Window for the TOA compensation
    'window_response_s': 2, # (s) Widow for response analysis
    'sound_speed': 1490, # m/s
    'whale_depth_m': 30, # Whale vocalizing depth (m)
    'whale_freq': 20, # Hz
    'whale_SL_dB': 189, # dB re 1uPa @1m
    'alpha_cable_fiber': 0.8, #
    'bulk_water': 2.29 * 10 ** 9,
    'H_multiple': 4, # Regime shift x water colum height ### For Svalbard replace by 8
}

# Define path to folder with data #TODO: download data  from XXXX
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

# Print the available datasets and the different fields
print(f'List of the datasets names: {list(all_dataset_configs.index)}')

# Select the dataset to work with 'OOISouthC1' or 'Svalbard' or 'MedSea'
dataset = 'OOISouthC1'
dataset_config = all_dataset_configs.loc[dataset]

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

# Part 1: Time compensate the labeled data
# 1) Load and condition the das data
print(f'Processing ... {os.path.basename(file)}')

if dataset_config['file_duration_s']-start_time_s < analysis_parameters['window_toa_compensation_s']:
    file_list = find_next_file(os.path.join(data_dir,file))
else:
    file_list = os.path.join(data_dir,file)

# Load the labeled file and the next to avoid edge effects in the time compensation.
tr, time, dist, fileBeginTimeUTC, metadata, selected_channels_m = load_data(dataset_config['interrogator'],
                                                                            os.path.join(data_dir,file),
                                                                            dataset_config['selected_channels_m'],
                                                                            sensitivity=dataset_config['sensitivity'])

fs_original, dx, nx, ns, gauge_length = metadata["fs"], metadata["dx"], metadata["nx"], metadata["ns"], \
    metadata["GL"]
print(f'File duration: {ns/fs_original}')


selected_channels = [int(selected_channels_m // metadata["dx"]) for selected_channels_m in
                     selected_channels_m]
# Adjust to ensure the range created by selected_channels is even length
if (selected_channels[1] - selected_channels[0]) % selected_channels[2] != 0:
    # Modify the first value so the length of the range becomes even
    selected_channels[0] -= 1

print(selected_channels)

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
print(f'Check dimensions: Dist {np.shape(dist)}, other {np.shape(position["depth"])}')

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

# Part 2: Response

# Get GL and LW
#GL = metadata['GL']
#LW =
LW = gauge_length

# 1) Calculate the distances and angles
# Depth difference z
z = np.abs(np.array(position['depth'])) - analysis_parameters['whale_depth_m']

# Horizontal distance to the closest channel hy
hy = whale_position['offset']

# Distance between apex and channels, hx
# Find the index of whale_apex_m in dist
ind_whale_apex = np.where(dist >= whale_position['apex'])[0][0]

# Correct the whale offset if 0, becomes depth-whale depth
whale_offset_m_temp = int(np.floor(np.array(position['depth'])[ind_whale_apex] - analysis_parameters['whale_depth_m']))
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

# 3)  Received level at the cable
# Propagation (geometrical spreading + interference from Lloyd's)
g_lme = 2 * np.sin(2 * np.pi * analysis_parameters['whale_freq'] * analysis_parameters['whale_depth_m'] * cos_phi / analysis_parameters['sound_speed'])
g_spherical =  1/r
TL_water_spherical_dB =  20 * np.log10(g_lme) + 20 * np.log10(g_spherical)

g_cylindrical = 1 / np.sqrt(r)
TL_water_cylindrical_dB =  20 * np.log10(g_cylindrical)

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
                                np.where(r_positive >= analysis_parameters['H_multiple'] * position['depth'][ind_whale_apex])[0][0]
    TL_water_regime_shift_dB[ind_regime_shift_positive:] = -10 * np.log10(
        r[ind_regime_shift_positive - 1]) + 20 * np.log10(
        g_lme[ind_regime_shift_positive - 1]) - 10 * np.log10(r[ind_regime_shift_positive:])
except:
    print('No regime shift')

# Conversion strain to Pa
G_strain_to_pa = 20 * np.log10(analysis_parameters['bulk_water'])

# Coupling between cable and fiber
g_coup_fiber_cable = (0.7 * np.cos(theta) ** 2 - 0.2 * analysis_parameters['alpha_cable_fiber'] * np.sin(theta) ** 2)

# Effect of the gauge length
g_GL = (np.sin((np.pi * analysis_parameters['whale_freq'] / analysis_parameters['sound_speed']) * gauge_length * np.cos(theta)) / (
            (np.pi * analysis_parameters['whale_freq'] / analysis_parameters['sound_speed']) * gauge_length * np.cos(theta))) * (
                   np.sin((np.pi * analysis_parameters['whale_freq'] / analysis_parameters['sound_speed']) * LW * np.cos(theta)) / (
                       (np.pi * analysis_parameters['whale_freq'] / analysis_parameters['sound_speed']) * LW * np.cos(theta)))

# Receive level on DAS (strain)
RL_DAS_spherical_dB = analysis_parameters['whale_SL_dB'] + TL_water_spherical_dB + 20 * np.log10(g_GL) + 20 * np.log10(
    g_coup_fiber_cable) - G_strain_to_pa

RL_DAS_cylindrical_dB =  analysis_parameters['whale_SL_dB'] + TL_water_cylindrical_dB  + 20 * np.log10(g_GL) + 20 * np.log10(
    g_coup_fiber_cable) - G_strain_to_pa

RL_DAS_regime_shift_dB = analysis_parameters['whale_SL_dB'] + TL_water_regime_shift_dB + 20 * np.log10(g_GL) + 20 * np.log10(
    g_coup_fiber_cable) - G_strain_to_pa

# DAS Response
DAS_response = 20 * np.log10(g_GL) + 20 * np.log10(g_coup_fiber_cable) #- attenuation_DAS * dist / 1000

# Receive level on DAS (dB Pa)
RL_DAS_Pa_dB = 20 * np.log10(np.mean(abs(trf_fk_align_call), axis=1) / 1e-6) - DAS_response + G_strain_to_pa

# Finding an offset to plot the theroretical curves
mean_resp = np.mean(abs(trf_fk_align_call), axis=1)
offset_curves = np.mean(20 * np.log10(mean_resp[mean_resp > 0] / 1e-6)) + 40

# Plot step by step  -------------------------------------
# Create a figure with 3 subplots
if dataset == 'OOISouthC1':
    width = 6
else:
    width = 5.45
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(width, 10), facecolor='none')

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

# Plot 2: Received level fiber
# Transmission loss and Lloyds Mirror effect
ax2.plot(dist / 1000, analysis_parameters['whale_SL_dB'] + TL_water_spherical_dB, 'steelblue', label=r'$RL_{dB \; re. 1\mu Pa} \; (r_T \to \infty)$')
ax2.plot(dist / 1000, analysis_parameters['whale_SL_dB'] + TL_water_cylindrical_dB - 30, 'steelblue', linestyle='dashed',
         label=r'$RL_{dB \; re. 1\mu Pa} \; (r_T \to 1) - 30 \;$ dB')
ax2.plot(dist / 1000, analysis_parameters['whale_SL_dB'] + TL_water_regime_shift_dB, 'k', label=r'$RL_{dB \; re. 1\mu Pa}}$')

# Labels
ax2.tick_params( labelsize=fontsize_vals)
ax2.set_ylabel('Rec. level (dB re. 1$\mu$Pa)', fontsize = fontsize)
ax2.yaxis.set_label_coords(-0.09, 0.5)

ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax2.set_xlabel('')  # Remove x-axis label for ax2

# Limits
ax2.set_ylim([95, 153])  # For the global plot
ax2.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax2.grid()

# Plot 2: DAS Response (strain)
# Coupling
ax3.plot(dist / 1000, 20 * np.log10(g_coup_fiber_cable), 'maroon', label='$G_{\kappa}$')
# Gauge length
ax3.plot(dist / 1000, 20 * np.log10(g_GL), 'indianred', linestyle='dashdot', label='$G_{GL}$')
# Combined
ax3.plot(dist / 1000, DAS_response, 'k', label=r'$H_{DAS_{\varepsilon}}$')

# Labels
ax3.set_ylabel(r'DAS response (dB re. 1$\mu\varepsilon$)', fontsize = fontsize)
ax3.yaxis.set_label_coords(-0.09, 0.5)
ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax3.set_xlabel('')  # Remove x-axis label for ax3
ax3.tick_params(labelsize=fontsize_vals)

# Limits
ax3.set_ylim([-25, 5])  # For the global plot
ax3.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax3.grid()

# Plot 3: Total RL reference in strain
#ax4.plot(dist / 1000, 20 * np.log10(np.mean(abs(noise), axis=1) / 1e-6), ### FOR OOI replace by dist_noise
#         label=r'$NL_{dB \; re. 1\mu \varepsilon} \;$ Measured', color='slategrey', alpha = 0.8, linewidth = 0.5)
ax4.plot(dist / 1000, 20 * np.log10(np.mean(abs(trf_fk_align_call), axis=1) / 1e-6),
         label=r'$RL_{dB \; re. 1\mu \varepsilon} \;$ Measured', color='orangered')

# Plot the regime shift
ax4.plot(dist / 1000, RL_DAS_regime_shift_dB, 'k', label=r'$RL_{dB \; re. 1\mu \varepsilon} \;$ Simulated, $S_{Pa \to \varepsilon} = 187.2$ dB')


# Labels
ax4.tick_params(labelsize=fontsize_vals)
ax4.set_ylabel(r'Rec. level (dB re. 1$\mu\varepsilon$)', fontsize = fontsize)
ax4.yaxis.set_label_coords(-0.09, 0.5)
ax4.set_xlabel('Distance (km)', fontsize = fontsize)

ax4.set_ylim([-105, -50])  # For the global plot
ax4.set_xlim([dist[0] / 1000, dist[-1] / 1000])
ax4.grid()

## FOR SVALBARD
#ax1.set_xlim([50 , dist[-1] / 1000])
#ax2.set_xlim([50, dist[-1] / 1000])
#ax3.set_xlim([50, dist[-1] / 1000])
#ax4.set_xlim([50, dist[-1] / 1000])
#ax2.legend(loc='upper right', fontsize = fontsize_legend)
#ax3.legend(loc='lower right', fontsize = fontsize_legend)
#ax4.legend(loc='upper right', fontsize = fontsize_legend)
#ax1.set_ylim([350, 0])
#manual_ticks = [50, 100, 150, 200, 250, 300]  # Adjust these values as needed

# FOR OOI South
#ax1.set_ylim([650, 0])
#manual_ticks = [100, 200, 300, 400, 500, 600]  # Adjust these values as needed

# FOR MedSEa
ax1.set_ylim([2500, 0])
manual_ticks = [500, 1000, 1500, 2000]  # Adjust these values as needed

ax1.set_yticks(manual_ticks)

# FOR ALL EXCEPT OOI --- Not all figs are the same if we do that
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