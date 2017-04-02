Scripts used to derive L and D functions from step-response measurements 
using the fluorescence polarization reporter (Fig. 3(C-D))
**************************************************************************
I performed two independent measurements of the step-up and step-down 
functions. The raw data and associated global fitting routines (see 
Computational Methods) are located in the following folders:
    stepDown_repeat_one_github/
    stepUp_repeat_one_B_github/
    stepUpDown_repeat_two_github/
Global fitting is performed on each step-up/step-down experiment separately
using the scripts in the corresponding directories. The two measurements
of L and two measurements of D are then merged into the four possible
combinations of L&D using the scripts in the merge_all_measured_L_and_D/ 
directory. The merged dataset is then used in simulations in Fig. 4 and 
Fig. 5.

This README file describes how to use the fitting pipelines to process
the measurements mentioned above. The fitting is done analogously in all 
folders (only the inputs and raw data files differ), so the steps are 
described only once. 

To generate L and D functions from polarization data data
*****************************
1. Run globalFit_steps_NonParamBoot.m with NUMBOOT=1 (no resampling of 
data, so do only one fit to all of the datapoints together) to generate
an output file containing the appropriate step-response functions and 
best-fit periods. 
Set TOSAVESTEPFUNS=1 to save the output into a timestamped .mat file. 
One such output file is provided in saved_data/ subfolder.
You can check goodness of fit by setting to TOTEST=1 in 
bestFitNonParBoot.m to plot fits on top of raw data. 
2. Do this for every dataset from which you want to get L and D funs (i.e.,
for every directory listed above).

To plot experimentally measured L and D functions (Fig. 3(C-D))
*************************************************
1. Open experimental_L_and_D_funs.m in merge_all_measured_step_funs.m.
2. Enter appropriate file names into the STEPUP_ONE_FILE, STEPDOWN_ONE_FILE
and STEPUPDOWN_FILE fields and corresponding directories. Note that two
directories contain a single step-up or step-down function measurement
(stepDown_repeat_one_github and stepUp_repeat_one_github), whereas 
stepUpDown_repeat_two contains measurements of both a step-up and a 
step-down function. Please enter the file names accordingly.
3. Run the file. This will generate the figures of L and D functions in
Fig. 3(C-D).

To generate the heatmap m as a function of L and D slopes (Fig. 4 - figSup3)
*********************************************************
1. Run mSlope_LD.m. This generates the heat map with the crosshair marking 
the best estimates of l and d based on my measurements.
2. The coordinates for crosshair marking l and d slope
intervals are computed by experimental_L_and_D_funs.m and are printed out
when that function is running. A copy of these values is saved in mSlope_LD.m.

To run and analyze simulations using L and D functions
******************************
Refer to step_fun_analysis_github/ in the parent directory to run 
simulations based on L and D measured from fluorescence polarization 
datasets, as described in Figure 4 and its supplements, Fig.3-figSup1, 
Figure 5 and Fig. 5-figure supplement. 


