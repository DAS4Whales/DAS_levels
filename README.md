# DAS Levels
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15008957.svg)](https://doi.org/10.5281/zenodo.15008957)
[![License](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-blue)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

<!-- SPHINX-START -->
## About DAS Levels 
This GitHub repository was created to support this paper [Link published paper].
We used this repository to investigate DAS received levels formalized as 
RL_(dB re. 1μϵ) = RL_(dB re. 1μPa) + S_(Pa → ϵ) + H_(DAS_ϵ)

where:
- RL_(dB re. 1μPa) is the received pressure level in dB re 1μPa,
- S_(Pa → ϵ) denotes the conversion term in dB between acoustic pressure and the strain measured by DAS, and
- H_(DAS_ϵ) denotes DAS's response in dB relative to strain.

Our work used labeled (apex, offset, start time) fin whale 20 Hz call dataset, available in open access through 
[Zenodo](https://doi.org/10.5281/zenodo.15008560)

This repository extensively uses functions from [DAS4Whales](https://github.com/DAS4Whales/DAS4Whales), a Python library to analyze 
Distributed Acoustic Sensing (DAS) data for marine bioacoustics.

Authors: 
- Léa Bouffaut, Ph.D.   
K. Lisa Yang Center for Conservation Bioacoustics, Cornell Lab of Ornithology, Cornell University. Ithaca (NY), USA

[Contact](mailto:lea.bouffaut@cornell.edu) | [Webpage](https://www.birds.cornell.edu/ccb/lea-bouffaut/)


```{note}
Please cite this code as:   
Léa Bouffaut (2025). DAS4Whales: A Python package to analyze Distributed Acoustic Sensing (DAS) data for marine bioacoustics (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.7760187
```

## Workflow 
To investigate the amplitude response of the DAS, the data undergoes multiple steps to account 
for instrument-specific variations in output, isolate fin whale 20~Hz calls and correct for 
the source-receiver geometry. 

### Labeled datasets
We recommend to download our openly accessible labeled dataset [openly accessible labeled dataset](https://doi.org/10.5281/zenodo.15008560),

Data were labeled using a custom annotation app: [DASSourceLocator](https://github.com/leabouffaut/DASSourceLocator). 
The annotator defined the apex (channel of the first time of arrival), offset (distance between the whale and the 
DAS at the apex), and first time of arrival based on the match between the curves. When the cable geometry was 
sufficiently asymmetric, the annotator could determine the source's side relative to the cable.

### Data processing
Data is pre-processed using the [DAS4Whales](https://github.com/DAS4Whales/DAS4Whales) Python package to 
harmonize temporal and spatial sampling. Each file is conditioned by:
1. loading and converting data to strain (interrogator-dependent),
2. resampling data in time and space and, 
3. applying a combined frequency-wave number ($f$--$k$) and bandpass filter
to isolate fin whale 20 Hz calls propagating in the water and sediments.
4. The data is time-compensated to measure DAS received levels using theoretical arrival 
times derived from the annotated whale position.
