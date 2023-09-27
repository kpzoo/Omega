# Angular reproduction numbers, Omega

Code to compute the angular reproduction number (Omega) and compare it to the effective reproduction number (R) and growth rate (r). The angular reproduction number provides estimates of transmissibility that are more robust than standard approaches to misspecified or time-varying generation time distributions.

Omega is introduced in the paper "Angular reproduction numbers improve estimates of transmissibility when disease generation times are misspecified or time-varying" by Kris V Parag, Benjamin Cowling and Ben C Lambert at [medRxiv 2022.10.19.22281255; doi: https://doi.org/10.1101/2022.10.19.22281255. This has been accepted at](https://royalsocietypublishing.org/doi/10.1098/rspb.2023.1664) in the Proceedings of the Royal Society B.

System Requirements

Matlab code should work with any standard Matlab installation and is tested on macOS v13.4.1 using Matlab v2022a. Dependencies are included in the main folder. R files tested on R version 4.2.2 using RStudio version 2023.03.1+446.

Instructions and installation

Run FigX.m to generate and reproduce the respective figure from the preprint. The main folder contains dependencies from other packages including EpiFilter, which is used to estimate reproduction numbers (see https://github.com/kpzoo/EpiFilter).

The omega R version folder contains scripts in R to compute Omega and R on custom datasets (omegaAliExample - which reproduces Fig 6 of the above paper which uses data from Ali et al (2020) Science 369(6507):1106-1109, doi: 10.1126/science.abc9004) or using simulated incidence curves (omegaSim and omegaSimVary). Outputs go the the results subfolder.

Both Matlab and R scripts are self contained so no external installations required. Run times of all scripts are of the order of minutes or faster.
