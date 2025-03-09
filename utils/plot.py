import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import numpy as np

def get_ticks_dist(selected_channels_m):

    # Round up min and round down max to nearest 5 km
    rounded_min = np.ceil(selected_channels_m[0]/1000 / 5) * 5
    rounded_max = np.floor(selected_channels_m[1]/1000 / 5) * 5

    # Generate ticks every 5 km
    ticks_dist = np.arange(rounded_min, rounded_max + 5, 5)

    return ticks_dist


def plot_save_TOA_compensated(trace, dist, window_toa_compensation_s,
                          file_begin_time_utc=datetime(2000, 1, 1, 1, 10, 10),
                          fig_size=(12, 10), v_min=None, v_max=None,
                              ticks_distance_km=None, Save=[True, 'img/', 'filename.png']):
    """
    Spatio-temporal representation (t-x plot) of the strain data

    Inputs:
    :param trace: a [channel x time sample] nparray containing the strain data in the spatio-temporal domain
    :param time: the corresponding time vector
    :param dist: the corresponding distance along the FO cable vector
    :param file_begin_time_utc: the time stamp of the represented file
    :param fig_size: Tuple of the figure dimensions. Default fig_size=(12, 10)
    :param v_min: sets the min nano strain amplitudes of the colorbar. Default v_min=0
    :param v_max: sets the max nano strain amplitudes of the colorbar, Default v_max=0.2

    Outputs:
    :return: a tx plot

    """
    fig = plt.figure(figsize=fig_size, num='TX')  # + file_begin_time_utc.strftime("%Y%m%d-%H%M%S"))
    shw = plt.imshow(abs(trace) * 10 ** 9, extent=[-window_toa_compensation_s, window_toa_compensation_s, dist[0] * 1e-3, dist[-1] * 1e-3, ], aspect='auto',
                     origin='lower', cmap='jet', vmin=v_min, vmax=v_max)
    plt.plot([0, 0], [dist[0] * 1e-3, dist[-1] * 1e-3], color = 'w')
    plt.plot([2, 2], [dist[0] * 1e-3, dist[-1] * 1e-3], color = 'w')

    plt.yticks(ticks_distance_km)
    plt.ylabel('Distance (km)')
    plt.xlabel('Time (s)')
    plt.xlim([-window_toa_compensation_s, window_toa_compensation_s])
    bar = fig.colorbar(shw, aspect=20)
    bar.set_label('Strain (x$10^{-9}$)')

    if isinstance(file_begin_time_utc, datetime):
        plt.title(file_begin_time_utc.strftime("%Y-%m-%d %H:%M:%S"), loc='right')

    # Save
    if Save[0]:
        print()
        plt.savefig(Save[1], bbox_inches='tight', transparent=True, dpi=180)
        plt.show(block=False)
        plt.clf()
        # plt.close('TX' + file_begin_time_utc.strftime("%Y%m%d-%H%M%S"))
        plt.close(fig)
    else:
        plt.show(block=False)
        plt.close(fig)