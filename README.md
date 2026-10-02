# SN-ORACLE

### *S*uper*N*ova *O*utcome *R*ates *A*cross *C*omprehensive *L*arge-scale single and binary stellar *E*volution models,

This population synthesis code is a post-processing tool that analyzes a large-scale input grid of detailed single and binary stellar evolution models to predict numbers and features of supernovae.

This repository contains the source code of SN-ORACLE. The necessary data tables are available at [Zenodo](10.5281/zenodo.22208968). 

## Science Aim

The last few years have seen an explosion (pun intended) in the number of observed supernovae: from less than 20 in the year 1980 to ~200 in 2000, to ~20,000 in 2020. We are thus living in a golden age of supernova observations (see [The Open Supernova Catalog](https://sne.space/) or the latest news on the [Rochester Academy of Science website](https://www.rochesterastronomy.org/supernova.html)), pushed forward by the operation of all-sky surveys, especially the Zwicky Transient Facility (Palomar Observatory) and the Legacy Survey of Space and Time (Vera C. Rubin Observatory).
Our goal is to bridge these observations with current state-of-the-art large-scale model grids of massive-star supernova progenitors. As most massive stars live in binaries, the parameter space of theoretical models is huge, and only recently were we able to produce large-scale comprehensive grids. 
With this code, we translate detailed stellar/binary evolution models into theoretical predictions on the relative rates and properties of the different supernovae produced. 
This is the first step to use the observed sample of supernovae to learn something new about their progenitor stars. 


## Requirements

This code requires Julia to run, and we suggest Jupyter Notebooks for running the simulations and analyzing the data.

### Data 

This code requires data from stellar/binary evolution models, typically produced with codes like MESA. 
The current version of this code relies on input data from Jin et al. ([2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...690A.135J/abstract), [2026](https://ui.adsabs.harvard.edu/abs/2026A%26A...707A..56J/abstract)) and Ercolino et al. ([2026a](https://ui.adsabs.harvard.edu/abs/2026A%26A...706A.169E/abstract)). Their data is, for now, partially hardcoded within the code. Future releases of SN-ORACLE will aim to make the code flexible and operational with user input using different grids. 
Feel free to contact me if you want to use your single and binary evolution grids, and I'll see if it can be done with little effort! 

### System requirements

Simulations and the following data analysis may deal with even GBs worth of data. Keep that in mind when executing this program.

## Installation 
- Clone/Download the repository to your working folder, say `user/SN-ORACLE/`.
- Download the data folder from Zenodo. You must have a directory `user/SN-ORACLE/Data/` in order for the code to run.
- `cd user/SN-ORACLE/` and then `julia install_SN_ORACLE.jl`: this will install all necessary packages and dependencies and pre-compile SN-ORACLE for use. It may take ~10 minutes.

## Use 
Open runner.ipynb with your favorite editor like VSCode. There, you will find a click-and-play example for running SN-ORACLE, reproducing the Fiducial Model in Ercolino et al. (2026a). A user's guide will be written in the future.

