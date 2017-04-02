Data for phase-resetting (PRC and wedge) and entrainment in Fig. 5 
and scripts used to perform a global fit to these datasets
************************************************************************

Data from experiments, plotted in Fig. 5 panels, is in the following files:
*********************

A. PRC: 2016-10-25_16.28.08_prcdata.mat
Data organized as follows (dp stands for dark pulse):
[dpDuration dpTime avg_phase_shift phase_shift_stderr phase_shift_stdev]

B. Wedge: 2016-10-25_23.24.46_wedgedata.mat
Data organized as follows (dp stands for dark pulse):
[dpDuration dpTime avg_phase_shift phase_shift_stderr phase_shift_stdev]

C. Entrainment: 2016-10-25_16.56.02_TeeTauData.mat
Data organized as follows  (tau = day length, T = driving period (22-26 hrs),
pkT = peak time after release into LL)
[tau T pkT pkT_stderr pkT_stdev]


Code used to generate global fit and all panels in Fig. 5
*********************************************************
fitMultipleDatasets.m