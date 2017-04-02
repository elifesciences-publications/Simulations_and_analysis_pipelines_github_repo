Scripts used for bifurcation analysis (fig. 3-figSup 1)
*************************************

To run simulations 
******************
1. Open runPhOsc_Boot_LinNLin_Freq.m and input appropriate parameters to 
drive the phase oscillator for specific driving periods, number of 
light-dark cycles, and day length (default=0.5*driving period).
2. Hit run. Once complete, the program will generate a time-stamped file
TIMESTAMP_freqSweep.mat that will contain your simulation results. The 
program will also save a timestamped version of the simulation script.

To plot simulation results
**************************
1. Open analyzePhOsc_Sims_Boot_LinNLin_freqSweep.m and input the simulation 
parametes (look them up in the timestamped simulation script if you forgot).
2. Hit run -- this will generate the bifurcation diagram.

To plot simulation results in Fig. 3-figSup1
********************************************
1. Run analyzePhOsc_Sims_Boot_LinNLin_freqSweep.m using 
TIMESTAMP_freqSweep_fine.mat as the input file. 