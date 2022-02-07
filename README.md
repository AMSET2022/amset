# Paper repository 

This repository contains the code to produce the all the simulation results and figures in our paper: **Empirical Bayes multistage testing for large-scale experiments**.

## R files

The folder ```R``` contains all R files. 

- ```epsest.func.R``` contains auxiliary function to estimate the proportion of non-null locations. 
- ```MSET.R``` implements fixed horizon comparison between our multistage empirical Bayes testing procedure and Optimizely's method using always valid p-values. 
- ```MSET_stagewise.R``` implements stagewise comparison between our multistage empirical Bayes testing procedure and Optimizely's method using always valid p-values. 
- ```AMSET.R``` implements adaptive multistage empirical Bayes testing procedure. 

## Data files

The folder ```Data``` contains simulation results and figures in our paper. 

- ```fixed_horizon_res``` contains simulation results for fixed horizon comparisons. 
- ```stagewise_res``` contains simulation results for stagewise comparisons. 

## Software installation
The package ```REBayes``` that provides Kiefer-Wolfowitz maximum likelihood estimation for mixture models requires the installation of the optimization package```Rmosek```. Details of installation can be found in Rmosek user guide (https://docs.mosek.com/latest/rmosek/install-interface.html).

 
