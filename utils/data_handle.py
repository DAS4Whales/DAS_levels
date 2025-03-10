import csv
import das4whales as dw
import numpy as np
import h5py
from datetime import datetime, timedelta
import os
from simpledas import simpleDASreader as sd
import pandas as pd

# Load annotations
def get_analysis_csv(annotation_csv):

    # Prepare lists
    files_list = []
    fs = []
    dx_oj = []
    dx = []
    dist_min = []
    dist_max =  []
    channel_min =[]
    channel_max = []
    channel_step = []
    gauge_length = []
    apex_list = []
    offset_list = []
    start_time_list =[]
    whale_side_list = []

    with open(annotation_csv, mode='r') as file:
        csv_file = csv.reader(file, delimiter=',')
        next(csv_file)  # Skip the header
        for row in csv_file:
            if not row:  # Check if the row is empty
                continue  # Skip empty rows
            files_list.append(row[1])
            fs.append(float(row[3]))
            dx_oj.append(float(row[4]))
            dx.append(float(row[5]))
            dist_min.append(float(row[6]))
            dist_max.append(float(row[7]))
            channel_min.append(int(row[8]))
            channel_max.append(int(row[9]))
            channel_step.append(int(row[10]))
            gauge_length.append(float(row[11]))

            apex_list.append(int(row[12]))  # Convert to integer
            offset_list.append(int(row[13]))  # Convert to integer
            start_time_list.append(float(row[14]))  # Convert to float
            whale_side_list.append(row[15])


    # Create the dictionary
    analysis_dict = {
        'TOA-compensated file': files_list,
        'fs export (Hz)': fs,
        'dx original (m)': dx_oj,
        'dx export (m)': dx,
        'dist min (m)': dist_min,
        'dist max (m)': dist_max,
        'channel min':channel_min,
        'channel max':channel_max,
        'channel step': channel_step,
        'gauge length (m)': gauge_length,
        'Whale apex (m)': apex_list,
        'Whale offset (m)': offset_list,
        'Start time (s)': start_time_list,
        'Whale side': whale_side_list
    }

    return analysis_dict

# Find files
def find_next_file(file_path):
    # Get the directory and filename
    directory = os.path.dirname(file_path)
    filename = os.path.basename(file_path)

    # List all files in the directory
    files = os.listdir(directory)

    # Sort the files
    sorted_files = sorted(files)

    # Find the index of the sought file
    index = sorted_files.index(filename)

    # Get the following file
    if index + 1 < len(sorted_files):
        file_path = [os.path.join(directory, filename),
                     os.path.join(directory, sorted_files[index + 1])]
    return file_path

# Load DAS data
def load_multiple_optasense_data(filename_list, selected_channels, metadata):
    """
    Load the DAS data corresponding to the input file name as strain according to the selected channels.

    Parameters
    ----------
    filename_list : str, list
        A string containing the full path to the data to load.
    selected_channels : list
        A list containing the selected channels.
    metadata : dict
        A dictionary filled with metadata (sampling frequency, channel spacing, scale factor...).

    Returns
    -------
    trace : np.ndarray
        A [channel x sample] nparray containing the strain data.
    tx : np.ndarray
        The corresponding time axis (s).
    dist : np.ndarray
        The corresponding distance along the FO cable axis (m).
    file_begin_time_utc : datetime.datetime
        The beginning time of the file, can be printed using file_begin_time_utc.strftime("%Y-%m-%d %H:%M:%S").
    """
    raw_data_time = []
    for filename in filename_list:
        if not os.path.exists(filename):
            raise FileNotFoundError(f'File {filename} not found')

        with h5py.File(filename, 'r') as fp:
            # Data matrix
            raw_data = fp['Acquisition/Raw[0]/RawData']

            # UTC Time vector for naming
            raw_data_time.append(fp['Acquisition']['Raw[0]']['RawDataTime'])

            # Check the orientation of the data compared to the metadata
            if raw_data.shape[0] == metadata["nx"]:
                # Data is in the correct orientation nx x ns
                pass
            elif raw_data.shape[1] == metadata["nx"]:
                # Data is transposed without loading in memory
                raw_data = raw_data[:,:].T

            # Now append the matrices to trace
            if raw_data.shape[1] == metadata["ns"]: # means only one file
                trace = raw_data
            else:
                trace = np.hstack((trace,raw_data))

    # Keep only the first time
    raw_data_time = raw_data_time[0]

    # Selection the traces corresponding to the desired channels
    # Loaded as float64, float 32 might be sufficient?

    trace = trace[selected_channels[0]:selected_channels[1]:selected_channels[2], :].astype(np.float64)
    trace = dw.data_handle.raw2strain(trace, metadata)


    # For future save
    file_begin_time_utc = datetime.utcfromtimestamp(raw_data_time[0] * 1e-6)

    # Store the following as the dimensions of our data block
    nnx = trace.shape[0]
    nns = trace.shape[1]

    # Define new time and distance axes
    tx = np.arange(nns) / metadata["fs"]
    dist = (np.arange(nnx) * selected_channels[2] + selected_channels[0]) * metadata["dx"]

    return trace, tx, dist, file_begin_time_utc

def load_multiple_asn_data(filename_list, selected_channels, metadata):
    dfdas = sd.load_DAS_files(filename_list, chIndex=None, samples=None, sensitivitySelect=-3,
                              userSensitivity={'sensitivity': metadata['scale_factor'],
                                               'sensitivityUnit': 'rad/(m*strain)'},
                              integrate=True, unwr=True)

    trace = dfdas.values.T
    trace = trace[selected_channels[0]:selected_channels[1]:selected_channels[2], :].astype(np.float64)

    # For future save
    file_begin_time_utc = dfdas.meta['time']
    
    # Store the following as the dimensions of our data block
    nnx = trace.shape[0]
    nns = trace.shape[1]

    # Define new time and distance axes
    tx = np.arange(nns) / metadata['fs']
    dist = (np.arange(nnx) * selected_channels[2] + selected_channels[0]) * metadata['dx']

    return trace, tx, dist, file_begin_time_utc

def load_data(interrogator, file, selected_channels_m, sensitivity=('default', None)):

    # Validate sensitivity input
    if not isinstance(sensitivity, tuple) or len(sensitivity) != 2:
        raise ValueError("Sensitivity must be a tuple of the form ('default', None) or ('custom', float).")

    mode, value = sensitivity

    if interrogator == 'optasense': # OOI data
        # 1) Get metadata
        if isinstance(file,list):
            metadata = dw.data_handle.get_acquisition_parameters(file[0], interrogator='optasense')
        else:
            metadata = dw.data_handle.get_acquisition_parameters(file, interrogator='optasense')

        # Select the channels
        selected_channels = select_channels_m_to_nb(selected_channels_m, metadata)

        # Loads the data using the pre-defined selected channels.
        if isinstance(file,list):
            tr, time, dist, fileBeginTimeUTC = load_multiple_optasense_data(file, selected_channels, metadata)
        else:
            tr, time, dist, fileBeginTimeUTC = dw.data_handle.load_das_data(file, selected_channels, metadata)

    elif interrogator == 'asn':
        # 1) Get metadata
        if isinstance(file,list):
            metadata = dw.data_handle.get_metadata_asn(file[0])
        else:
            metadata = dw.data_handle.get_metadata_asn(file)

        # 2) Adjust the sensitivity as needed in the metadata so it can be custom-set
        if mode == 'default':
            if value is not None:
                raise ValueError("For 'default', the second element of the tuple must be None.")

        elif mode == 'custom':
            if not isinstance(value, (float, int)):
                raise ValueError("For 'custom', the second element of the tuple must be a float or int.")
            metadata['scale_factor'] = value

        else:
            raise ValueError("The first element of the tuple must be 'default' or 'custom'.")

        # 3) Select the channels
        selected_channels = select_channels_m_to_nb(selected_channels_m, metadata)

        # 4) Load the data using the pre-defined selected channels.
        tr, time, dist, fileBeginTimeUTC = load_multiple_asn_data(file, selected_channels, metadata)

    return tr, time, dist, fileBeginTimeUTC, metadata, selected_channels_m

def select_channels_m_to_nb(selected_channels_m,metadata):
    selected_channels_m[2] = round(selected_channels_m[2] / metadata["dx"]) * metadata["dx"]
    selected_channels = [int(selected_channels_m // metadata["dx"]) for selected_channels_m in
                         selected_channels_m]

    # Adjust to ensure the range created by selected_channels is even length
    if (selected_channels[1] - selected_channels[0]) % selected_channels[2] != 0:
        # Modify the first value so the length of the range becomes even
        selected_channels[0] -= 1
    return selected_channels







