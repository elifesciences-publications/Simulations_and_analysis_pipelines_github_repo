2017-June-27, Eugene Leypunskiy (leypunskiy@uchicago.edu)
******************************

This directory contains scripts and functions used to analyze data and perform simulations described in "The cyanobacterial circadian clock tracks midday in vivo and in vitro" by Leypunskiy et al.

Below is a map describing where you can find scripts for specific figures/analyses. Each directory contains a readme.txt file that describes how the files are organized and which files perform which analyses.

-------------------------------------------------------------------------

Usage Notes
***********
1. A lot of the scripts rely on accessory files provided in dependencies_github/ folder. Please add those files and subfolders to your MATLAB path before running any scripts. 

2. A number of scripts contain variables that control whether or not a given function saves output, displays progress messages, makes plots or exports figures as PDF files. These variables have names such as TOPLOT, TOEXP, TODISP, TOSAVE, etc. Set their values to 1 to perform the relevant action (e.g., TOPLOT=1 turns plotting on).

3. On some systems, figure export fails because of an error in export_fig.m*. If you get an error message, try turning off figure export in the relevant file (e.g., by setting TOEXP=1) and re-running the code. Generally, errors in figure export should not affect the rest of the functionalities: the code should run, save output results appropriately and generate figures (but not export them). You can manually save the figures by using Matlab's dialog boxes.

*export_fig.m is a useful figure exporting package by Oliver Woodford, used and distributed according to the BSD license.

-------------------------------------------------------------------------

Figure 2 (including supplements)
********
Refer to in_vitro_entrainment_github/. 

-------------------------------------------------------------------------

Figure 3C-D
********
Refer to step_fun_from_polarization_github/.

Figure 3E-F
********
Refer to step_fun_pkaic_github/ for analysis of step-response functions measured by monitoring KaiC phosphorylation using SDS-PAGE.

Figure 3 - figure supplement 1
******************************
Refer to step_fun_analysis_github/bifurcation analysis/.

-------------------------------------------------------------------------

Figure 4 (incl. figure supplements 1A-B, 2A, 2C, 3)
********
Refer to entrainment simulations in step_fun_analysis_github/entrainment simulations/.

Figure 4 - figure supplement 2B
*******************************
Refer to step_fun_pkaic_github/ for analysis of step-response functions measured by monitoring KaiC phosphorylation using SDS-PAGE.

Figure 4 - figure supplement 3
******************************
Refer to step_fun_from_polarization_github/readme.txt. 
Run mSlope_LD.m

-------------------------------------------------------------------------

Figure 5 
********
Refer to in_vivo_phase_oscilllator_model_fit_github/.

Figure 5 - figure supplement 
****************************
Refer to step_fun_analysis_github/PRC simulations/.

-------------------------------------------------------------------------

Figure 6 (incl. figure supplements 1-2)
********
Refer to geometric_resetting_github/.

Figure 6 - figure supplement 3
******************************
Refer to dawn_dusk_genes_mVals_github/.

-------------------------------------------------------------------------

Mathematical Appendix Figure (illustrations in A-D)
****************************
Refer to step_fun_pkaic_github/.