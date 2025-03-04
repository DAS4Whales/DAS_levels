Svalbard DAS data 

DAS data -- available in data/
This data was recorded on 20220822 on the Svalbard DAS fiber optic cable out of Longyearbyen, Norway in the North Atlantic Ocean(https://doi.org/10.3389/fmars.2023.1130898;  https://doi.org/10.18710/Q8OSON). It was instrumented with an Alcatel OptoDAS that covers 137.9 km.

Project PI: Martin Landrø, CGF - NTNU
Gauge length: 8 m
Sampling frequency: 625 Hz
Channel spacing: 4 m
Sediment types: Sludge (see https://www.ngu.no/en/geologiske-kart)
Date: 20220822-123957UTC

Cable positions -- available in DASPositions/ 
Picked from Nautical charts (https://norgeskart.no/ base map layer: Nautical charts)
Format: channel / Lat / Lon / depth (m)

Manual analysis -- available in annotations/
Conducted in September 2024 by Léa Bouffaut

This analysis was conducted using a custom-made app (DASSourceLocator) that enables channel selection, shows the outputs of fin whale 20 Hz cross-correlation between a file and a synthetic call. It displays peaks above a pre-defined threshold on a t-x plot. The user can then manually adjust for each analyzed file
- whale offset (m): 1500 m
- whale apex (m): 79.0 km
- begin time (s): 3.8 s

to match the corresponding theoretical TOAs (calculated for both left and right whale positions) to the data. Whale depth is fixed at 30 m and c=1490 m/s for all analyses. Best match is saved into a recap annotation table.

Here: 
- Data loaded using:
	load_svalbard_das_data(file, selected_channels, metadata)
- selected_channels_m = [40000, 100000, 8] m
- Fk filter: cs_min=1350, cp_min=1450, cp_max=3300, cs_max=3450,
- Synthetic signal: fmin=14, fmax=24, duration=0.70
- X-corr threshold = 0.5 (x-corr normalized by max of the autocorr of the synthetic signal)
