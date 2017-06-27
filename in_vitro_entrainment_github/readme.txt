Scripts used to process data from in vitro entrainment experiments
******************************************************************

This directory contains scripts used to perform fits to SDS-PAGE and 
polarization datasets in Figure and its supplements and to make
the corresponding figures. 

Figure 2C
********
fluorescence_polarization/fit_PKaiC_Polarization_diffPeriod_May2017.m
Performs a global fit of entrainment datasets collected using the
fluorescence polarization probe and overlays the resulting dayTime vs. day
length curves on top of SDS-PAGE data and in vivo data (as in Fig. 2C).

Figure 2 - figuSup 1
********************
pkaic_vs_polarization_github/fit_PKaiC_Polarization_joint.m
performs a global fit to the polarization and SDS-PAGE datasets in 
Fig2-figSup1 and generates the figure.