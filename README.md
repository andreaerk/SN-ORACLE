# SN-ORACLE
## SuperNova Outcome Rates Across Comprehensive Large-scale single and binary stellar Evolution models,

This population synthesis code is a post-processing tool which analyzes a large scale input grid of detailed single and binary stellar evolution models to predict numbers and features of supernove.

This repository contains the source code of SN-ORACLE. The necessary data tables are available at Zenodo. 

## Science Aim
The last few years have seen an explosion (pun intended) in the number of observed supernovae. From <20 in 1980 to ~200 in 2000, to ~20,000 in 2020. We want to bridge these useful observations with current state-of-the-art large-scale model grids of massive-star supernova progenitors. As most massive stars live in binaries, the parameter space of theoretical models is huge, and only recently we were able to produce large-scale comprehensive grids! 
With this code, we translate detailed stellar/binary evolution models into theoretical predictions on the relative rates and properties of the different supernovae produced. 
This is the first step to use the observed sample of supernovae to learn something new about their progenitor stars. 


## Requirements
### Data 
This code requires data from stellar/binary evolution models, typically produced with codes like MESA. 
Here, input data is already provided from Jin et al. (2024, 2026) and Ercolino et al. (2026). 
### Execution 
This code requires Julia to run, and we suggest Jupyter for running the simulations and analyzing the data. 
### System requirements
Simulations and the following data analysis may deal with even GBs worth of data. Keep that in mind when executing this program.

## Installation 
- Clone the repository
- julia install_SN_ORACLE.jl 
then, use runner.ipynb for examples to use the code.

