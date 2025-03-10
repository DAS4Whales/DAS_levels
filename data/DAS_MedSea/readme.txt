Med Sea DAS data 

DAS data -- available in data/
This data was recorded on 20230922 on the MEUST fiber optic cable out of Toulon, France in the Mediterranean Sea (https://doi.org/10.1121/10.0004129). It was instrumented with an Alcatel OptoDAS that covers 52 km.

Project PI: Anthony Sladen @Géoazur
Gauge length: 4 m
Sampling frequency: 72.7 Hz
Channel spacing: 4 m
Sediment types: Sand/Mud (see https://data.shom.fr)
Date: 20230922-091304UTC

Cable positions -- available in DASPositions/
Shared by project PI (picked from nautical charts). Channel positions are corrected for the best correspondance at locations that can be identified in the data: entry in the water and sharp changes of directions (the W between 40-50 km).
Format: channel / Lat / Lon / depth (m)


Manual analysis -- available in annotations/
Conducted in September 2024 by Léa Bouffaut

This analysis was conducted using a custom-made app (DASSourceLocator) that enables channel selection, shows the outputs of fin whale 20 Hz cross-correlation between a file and a synthetic call. It displays peaks above a pre-defined threshold on a t-x plot. The user can then manually adjust for each analyzed file
- whale offset (m): 2375
- whale apex (m): 35.35 km
- begin time (s): 0 s

to match the corresponding theoretical TOAs (calculated for both left and right whale positions) to the data. Whale depth is fixed at 30 m and c=1490 m/s for all analyses. Best match is saved into a recap annotation table.

How is the data processed for the manual analysis:
- ASN data loaded using: 
	tr, fileBeginTimeUTC, metadata = fct.load_ASN_DAS_file(list_file[0]) #ASN
- selected_channels_m = [3050, 53050, 8] m
- Fk filter: cs_min=1350, cp_min=1450, cp_max=3300, cs_max=3450,
- Synthetic signal: fmin=14, fmax=24, duration=0.70
- X-corr threshold = 0.5 (x-corr normalized by max of the autocorr of the synthetic signal)


