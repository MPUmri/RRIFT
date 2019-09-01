# Reference region and input function tail (RRIFT) method

## Introduction

This repository contains code for the manuscript:

Ahmed, Z., & Levesque, I. R. (2019). [Pharmacokinetic modeling of dynamic contrast‐enhanced MRI using a reference region and input function tail.](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27913) Magnetic Resonance in Medicine

A preprint is available in [this gitlab repository](https://gitlab.com/notZaki/rrift).

The manuscript uses the constrained extended reference region model (CERRM) which is implemented in `./mfiles/CERRM.m` while the reference region and input function tail (RRIFT) equation is implemented in `./mfiles/RRIFT.m`. 
The linearized extended Tofts model (ETM), which is used for comparison, is in `./mfiles/Tofts_LLSQ.m`.

A short example on simulated data is in `a00_quickExample.m`.

## Simulations

The simulation scripts have filenames starting from `b00` to `b05`. 
Since some of the simulation steps are time consuming (~30 minutes), the resulting `.mat` files have been uploaded into [an OSF repository](https://osf.io/wr3kf/files/). 
If the four .mat files (begining with the filename `sim`) are downloaded from the OSF repository and placed in the `data` folder, then all of the even-numbered scripts can be skipped---i.e. you can run `b01`, `b03`, and `b05` to produce the figures in the manuscript.

### Overview of simulation scripts

- `b00_makeSimMap.m` makes the 2D virtual phantom which contains 100 parameter combinations. This is used as a starting point for all the simulation.
    + The output file `simMap.mat` is in the OSF repository and should be placed in the `data` directory
- `b01_sketchOverview.m` will produce **Fig. 1** in the manuscript
    + This script also performs a quick comparison with the reference tissue plus vessel (RTPV) technique and displays the results in the console/terminal
- `b02_mainSimAnalysis.m` will fit RRIFT and the ETM to the simulated data under a range of noise levels and temporal resolutions
    + The output file `simResults.mat` is in the OSF repository and should be placed in the `data` directory
- `b03_mainSimFigures.m` will produce **Figs. 2 & 3** in the manuscript
- `b04_secondarySimAnalysis.m` is similar to `b02` except the reference tissue parameters are also allowed to vary. For speed, only a temporal resolution of 15 s and noise of 0.02 mM is used, but this can be changed in the code
    + The output files are: `simResultsTRes15-varKtRR.mat` and `simResultsTRes15-varVeRR.mat`. They should be placed in the `data` directory.
- `b05_secondarySimFigures.m` will produce **Fig. 4** in the manuscript

## In-vivo evaluation

### Obtaining the data

The data can be downloaded from the cancer imaging archive under the [TCGA-GBM collection](https://wiki.cancerimagingarchive.net/display/Public/TCGA-GBM), but since only a small subset of collection contains DCE-MRI data, the subset has been [uploaded online on the OSF](https://osf.io/wr3kf/).

The OSF repository can be used in one of two ways:

1. Download the DICOM data in the `TCGA-GBM` folder from [the OSF repository](https://osf.io/wr3kf/files/), and then use the script `x01_dicomReader.m` to perform the T1 mapping and convert signal to concentration.
2. Alternatively, the results from the previous step can be downloaded by selecting the `TCGA-GBM-MAT` folder from [the same OSF repository](https://osf.io/wr3kf/files/) and clicking on the `Download as zip` button. 
  Unzip the contents of the downloaded `.zip` into the `./data/TCGA-GBM-Mat` folder---make the folder if it doesn't exist.

After either of the above two options, the final directory tree should resemble:
```
RRIFT/
├── data/
│   ├── Results/
│   ├── TCGA-GBM-Masks/
│   └── TCGA-GBM-Mat/
│       ├── Ct/
│       ├── DCE/
│       ├── hdr/
│       └── T1/
├── mfiles/
:
```
where the subfolders contain:  
- `Ct/` - the concentration-time data
- `DCE/` - the signal-time data (i.e. acquired DCE-MRI data)
- `hdr/` - the DICOM headers for the DCE-MRI data
- `T1/` - the T1 map computed from the variable flip angle data

### Overview of in-vivo evaluation code

For convenience, the output of most of the scripts is included within the github repository. 
This means that the `c01` script can be skipped, and all the other scripts should work without any additional steps.

- `c01_preprocessDCE.m` grabs the relevant information---e.g. AIF, tumour curves, muscle curve---and saves them as a `.mat` file for future steps
    + The output is stored in `./data/TCGA-GBM-Results/c01_preprocessed/`
- `c02_doRRIFT.m` fits RRIFT and the ETM to the in-vivo data and saves the results
    + The output is stored in `./data/TCGA-GBM-Results/c02_postprocessed/`
- `c03_showResults` makes **Figs. 5, 8, S1, & S2** while additional features, such as noise, are reported in the console 
- `c04_showResults_SinglePatient` was used for **Figs. 6 & 7**
- `c05_doRRIFT_smallerTails` looks at the effect of reducing the tail duration, and produces a figure which was ultimately not included in the manuscript
- `c06_doRRIFT_downsampled` downsampled the in-vivo data and fits RRIFT and the ETM
    + The output is stored in `./data/TCGA-GBM-Results/c06_downsampled/`
- `c07_showResults_downsampled` makes **Fig. 10**
- `c08_showMaps_downsampled.m` makes **Fig. 9** 
