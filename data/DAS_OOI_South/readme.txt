OOI South DAS data 

DAS data -- available in data/
This data was recorded on 20211102 on the OOI South fiber optic cable out of Pacific City (OR), USA in the North East Pacific Ocean (https://doi.org/10.1121/10.0017104). It was instrumented with an OptaSense QuantX that covers 97 km. The complete dataset is available in open access (https://doi.org/10.58046/5J60-FJ89).

Project PI: William Wilcock, University of Washington
Gauge length: 51 m
Sampling frequency: 200 Hz
Channel spacing: 2 m
Sediment types: Sand/Mud (see https://data.shom.fr)
Date: 20211102-190014UTC

Cable positions -- available in DASPositions/ 
Shared by OOI in the online repository.
Format: channel / Lat / Lon / depth (m)

Manual analysis -- available in annotations/
Conducted in September 2024 by Léa Bouffaut

This analysis was conducted using a custom-made app (DASSourceLocator) that enables channel selection, shows the outputs of fin whale 20 Hz cross-correlation between a file and a synthetic call. It displays peaks above a pre-defined threshold on a t-x plot. The user can then manually adjust for each analyzed file
- whale offset (m): 340 m
- whale apex (m): 31.95 km
- begin time (s): 18 s

to match the corresponding theoretical TOAs (calculated for both left and right whale positions) to the data. Whale depth is fixed at 30 m and c=1490 m/s for all analyses. Best match is saved into a recap annotation table.

Here: 
- Data loaded using the DAS4Whales functions
- selected_channels_m = [15000, 65000, 8] m
- Fk filter: cs_min=1350, cp_min=1450, cp_max=3300, cs_max=3450,
- Synthetic signal: fmin=14, fmax=24, duration=0.70
- X-corr threshold = 0.2 (x-corr normalized by max of the autocorr of the synthetic signal)
- a bandpass filter is added to the NorthC2 data: tr = dw.dsp.bp_filt(tr, fs, 14, 30)
