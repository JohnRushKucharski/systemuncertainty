# System Uncertainty

This repository provide the code and data needed to reproduce the _Uncertainty decomposition to understand the influence of water systems model error in climate vulnerability assessments_ paper results.

## Steps to reproduce the papers results.
1. Run ‘error_model/delta_outflow_error_model.R’
2. Run ‘error_model/‘delta_pumping_error_model.R’  

Steps 1-2 produce output files in the ‘/error_model/’, ‘/error_model/simulated.files.outflows/‘, and ‘/error_model/simulated.files.pumping/‘ directories.

3. Run ‘zip_ensemble_simulations.py’

The ’zip_ensemble_simulations.py’ file and an intermediate file that reformats the the error model outputs (from step 1&2), preparing them for the sobel analysis (in step 4) it outputs the reformatted files in the '/sobol_analysis/ensemble-simulations-flood/' and '/sobol_analysis/ensemble-simulations-pumping/' directories.

4. Run ‘sobel_analysis/main_sobel.py’

Step 4 produces (and further manipulates) output files in the '/sobol_analysis' directory.

#### Optional Steps:

5. Run ‘error_model/Figure2_code.R’ file to reproduce figure 2

The ‘error_model/Figure2_code.R’ file writes its output to the ‘error_model/figs/‘ diretory.

**NOTE**: all file paths accessed in the R and Python scripts are relative to directory containing the script being run. For the scripts to run the file directory structure must be preserved. 

In a few cases these path variables may need to be adjusted (based on the configuration of the user’s system). For instance, if the default working directory for R is not the folder containing the .R file the first path assignment in each .R file will need to be adjusted. In case this becomes necessary, a list of path assignment locations in each scripts is provided below. 

delta_outflow_error_model.R path assignment on lines: [5, 6, 149, 177]
delta_pumping_error_model.R path assignment on lines: [5, 6, 130, 155]
zip_ensemble)simulations.py path assignment on lines: [11, 12, 23, 26]
main_sobol.py path assignments on lines: [10, 50, 61, 65, 70, 74, 77, 99]
(Note: if the file directory structure is maintained it is unlikely that the later path assignments in these file will need to be adjusted).