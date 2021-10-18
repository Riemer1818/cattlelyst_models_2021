## This is a piece of written text detailing how to navigate this gitlab folder and the meaning of the files. 

# Final_scripts: 
                all thesis functions, scripts and datafiles can be found purpose of the study ordered by folder name. 
# ODE_solver_errors: 
                all files needed to study negative and positive simulations, and find the relationship between parameters and these simulations,     also to be found in Negative_concentrations\Final_scripts.
#  Powerpoints: 
                Weekly updates from october on in powerpoint format, detailing progress.


# Final_scripts (DETAILED)
 ## Figures\Final_scripts
**What**: Folder names detail the figure to reproduce
**Organisation**: Folders are organized in a MATLAB-proof form, such that when the whole folder is downloaded, the figure can be computed. 
### Files\Figures:
Figure_XXXX correspond to the scripts to be run to compute the figure. 
Data# are datafiles for the different strains tested: 

        Data1: data for P. stutzeri YZN-001
        Data2: data for P. stutzeri XL-2
        Data_am_nia_ZN1: data for ammonia and nitrate strain P. putida ZN1
        Data_am_nii_ZN1: data for ammonia and nitrite strain P. putida ZN1
        Data_am_nia_SDU10: data for ammonia and nitrate strain P. stutzeri SDU10
        Data_am_nii_SDU10: data for ammonia and nitrite strain P. stutzeri SDU10
Output_XXX are files that hold the simulation data to be plotted
time_Volume_growth_matrix_##: contains the information to simulate the model with fixed volume.
 
        11 corresponds to P. stutzeri YZN-001
        22 corresponds to P. stutzeri XL-2
#### Determine_kg_and_figs\Figures:
Slightly different structure; contains Script_to_determine_kg, if this script is run, relationship between kg and ammonia steady-state   is       determined, result is plotted. 

## Generate_LHS_parameters\Final_scripts
**What**: For optimization procedure random parameter values were sampled by LHS. 
**Organisation**: create_parameters calls lhsdesign_modified to sample N parameter values. 
 
## Optimization_fixed_volume\Final_scripts
**What**: Files to run LHS and Simulated Annealing optimization procedure.
**Organisation**: in a MATLAB-proof form, datafiles, parameter sets, outputs, functions, scripts
### datafiles
Data1, Data2, time_Volume_growth_matrix_11, time_Volume_growth_matrix_22 (see line 16-27)
### parameter sets
        Folder: p_sets_SA_gas: 25 parameter sets resulting from simulated annealing. 
        parameters_27_11: 50,000 parameter sets used to optimize the model. 
        good_p_sets_gas_###: parameter sets corresponding to ## best fitting simulations. 
### Outputs: contain simulation data
        Folder: Output_50000_p_sets contains 5 files that contain simulation data (score, score_ij,sim_all) for 10,000 parameter sets each. 
        Output_50000_10_12: simulation data merged from folder Output_50000_p_sets
### Functions: run the optimization
        run_HNAD_server: used to run the LHS optimization, calls  optimise_HNAD_sim
            optimise_HNAD_sim is the optimization function and calls the model with fixed volume: HNAD_sim_Vol_optimize
        run_SA_gas: used to run the simulated annealing optimization, calls HNAD_sim_SA
            HNAD_sim_SA is the optimization function and calls the model with fixed volume: HNAD_sim_Vol_optimize
        HNAD_sim_Vol_optimize is the model function with fixed volume
### script
        order_simulations: script to order simulations according to score and isolate parameter values generating good simulations.

## Volume\Final_scripts
**What**: files to optimize volume parameters separately
**Organisation**: (MATLAB-proof) Separate folder detailing the optimization of ammonia, nitrate-transport and mumax parameters
### datafiles
        Vol_data11: population volume data corresponding to Data1
        Vol_data22: population volume data corresponding to Data2
        ammonia_uptake: change [ammonia] over time according to Data1
        nitrate_uptake: change [nitrate] over time according to Data2
        data1, data2
### parameter sets
        p_sets_vol_3.6_1_1: final mumax and KN values for growth YZN-001 (log-scale)
        p_sets_vol_3.6_2_2: final mumax and KN values for growth XL-2 (log-scale)
### functions
        run_SA_volume: runs simulated annealing to establish Monod parameter values, calls vol_opt (takes 1 min to run)
        vol_opt: function that calls the respective model function to simulate Monod, and computes the sum_squared_error
        monod1: function to simulate population volume over time for YZN-001
        monod2: function to simulate population volume over time for Xl-2
### SA_re_optimize_transport_and_mu
        folder that holds the files to run the simulated annealing to re-optimize parameters for transport and mumax. 
        run_SA_vol: runs simulated annealing with initial parameter estimate (p_sets_SA_vol_initial_mumax1), calls HNAD_SA_vol
        HNAD_sim_vol_SA: computes the score based on the provided parameter estimates and calls model function: HNAD_SA_vol
        HNAD_SA_vol: model function with experiment ID switch to simulate growth. 

## Negative_concentrations\Final_scripts
**What**: folder holds files to study negative concentrations 
**organisation**: Matlab-proof, files to determine Q, simulate negative concentrations, output files, data, and parameters.
### datafiles and Output
data1, data2, time_Volume_growth_matrix_11, time_Volume_growth_matrix_22 (see lines 16-27 for details)

        Q_values_neg_simulations_1: Q values corresp. to negative simulations.
        Q_values_pos_simulations_1: Q values corresp. to positive simulations. 
        Output_negative_7_98_best_1: simulation data for negative simulations. 
        Output_positive_91_98_best_1: simulation data for positive simulations. 
### parameter values
        p_sets_negative_7_98_best: from the 98 best parameter sets, the 7 parameter sets causing negative simulations
        p_sets_positive_91_98_best: from the 98 best parameter sets, the 91 parameter sets causing positive simulations
### functions, can be used to simulate positive simulations and negative simulations
        run_HNAD_model runs the model with provided parameters and saves the output, calls simulate_HNAD_sim
        simulate_HNAD_sim runs the model as frequently as the number of parameter sets provided, calls HNAD_sim_Vol 
        HNAD_sim_Vol model file with fixed volume. 
### scripts, can be used to determine Q
        determine_Q, file that determines Q values for multiple parameters and can be used to save the output. Corresponding derivation can be found in thesis report. 

## Mean_simulation\Final_scripts
**What**: folder contains files to obtain stacked simulations
**organisation**: Matlab-proof, datafiles, parameter sets, functions, simulation output. 
### datafiles and output files
data1, data2, time_Volume_growth_matrix_11,time_Volume_growth_matrix_22 (see lines 16-27 for details)

        Simulation_##_best_parameters: stacked simulation output with good_p_sets_gas_##
### parameter values
        good_p_sets_gas_## (see line 45 for more information)
### functions
        mean_HNAD_sim: simulates the system with provided parameters multiple times and stacks+concatenates simulations. 
        HNAD_sim_vol_optimize: model file with fixed volume. 
### script
        Mean_stdv_simulation: calls mean_HNAD_sim and saves the stacked simulation. (max I = 100) 

## sensitivity_analysis_full_model\Final_scripts
**What**: folder that contains files to conduct the full_model sensitivity analysis + generate figure
**organisation**: Matlab-proof, best parameter set, data files, 1 script, 2 functions
### datafiles
        data1, data2 (see lines 16-27 for details)
### parameter set  
        p_sets_best: final best parameter set. 
### functions
        sensitivity HNAD: function computes score and calls the model function, HNAD_sim_Vol
        HNAD_sim_Vol: model function with experiment ID to switch between growth functions and option for kg
### script
        script that runs the model with perturbed parameter sets to eventually plot relative sensitivity coefficients for 1%, 3%, 5%, and 10% parameter changes. 

## Limit_N2O_production\Final_scripts
**What**: OAT sensitivity analysis and pair-wise parameter perturbations to limit N2O Limit_N2O_production
**Organisation**: two different folders that contain the analyses, matlab-proof
### 3D mesh plots and percentages\Limit_N2O_production
#### datafiles and outputs
data1, data2 (see lines 16-27), these datafiles are redundant for the analysis.

        Output_N2O_SA_HNAD_sim_vol_#_#_1 (# (# can be 4,6,8,9,22,27 and correspond to parameters: Vmax4a, Vmax5a, Vmax6, Vmax7, Tmax2a, Tmax3b.) 
        data3, data for another P. stutzeri KT2, not used in the analysis, left in here to explain why the Output files have: score__NH2OH and sim_NH2OH, these files were not used, so please ignore. 
#### parameter values
        parameters_N2O_sensitivity_log_#_#, (# see line 147) 
#### .figs
        3D_mesh_#_#, details for # see line 147
#### functions
        create_parameters_N2O, by means of LHS design generate new parameter combinations of p_sets_best with two of parameters 4,6,8,9,22,27 altered. 
        run_HNAD_N2O, function that is supplied with changeing parameter pair and calls optimise_HNAD_sim_N2O
        optimise_HNAD_sim_N2O,performs model simulation and scores it, calls HNAD_sim_Vol
        HNAD_sim_Vol is model function with kg and experiment ID-wise growth switch
#### scripts
        mesh_plots_N2O_production, file that generates mesh plots, and percentages to compare mesh plots. 
### OAT_sensitivity_analysis\Limit_N2O_production
#### parameter values
        p_sets_best
#### .figs
        Figure_parameter_sensitivity_N2O_production_##, specify strain: XL-2 or YZN-001
#### functions
        sensitivity N2O, function that simulates the model and scores it accordingly, calls model function: HNAD_sim_Vol
        HNAD_sim_Vol, similar as described in line 160
#### script
        OAT_sensitivity_analysis_N2O_production script that calls sensitivity_N2O and provides the model parameter sets with 1 perturbed parameter. it consequently plots the absolute effect of the perturbation on the y-axis. 

## Simulate_heaviside_function\Final_scripts
**What**: contains files to run the model with the heaviside function
**organisation**: datafiles, 1 output file, best parameter set, functions
### datafiles and output
data1, data2 (see lines 16-27)

        Output_best_heaviside_1, output with best heaviside function
### parameter values
        p_sets_best final model parameter set
### functions
        run_HNAD_NO2, function to call optimise_HNAD_sim_NO2, and save output
        optimise_HNAD_sim_NO2, function that simulates model: run_HNAD_NO2, and scores the output
        HNAD_sim_Vol_heaviside: model function with heaviside function and kg. 

## Simulate_mixed_conditions\Final_scripts
**What**: folder that contains all files to simulate mixed conditions with functions adapted according to ZN1, and SDU10 data.
**organisation**: best parameter set, 5 functions
### parameter values
        p_sets_best, final model parameter set
### functions
        run_HNAD_mix, function that is used to simulate mixed conditions, specify which condition should be approximated here. 
        sim_HNAD_mix_#strain#, specify strain: SDU10 or ZN1, function to run the simulation and call model functions
        HNAD_sim_mix_#strain#, specify strain: SDU10 or ZN1, model function adapted for the specific conditions, i.e. nitrogen source preference

## Growth_different_ammonia0\Final_scripts
**What**: folder that contains all files to simulate conditions with varying initial ammonia concentration. 
**organisation**: parameter set, functions
### parameter values
        p_sets_best, final model parameter set
### functions
        run_HNAD_ammonia, used to simulate the model and save the output, calls sim_HNAD_am
        sim_HNAD_am, contains the initial conditions and is used to simulate model: HNAD_sim_am
        HNAD_sim_am, model file with kg and nitrogen switches. 

_If you have any further questions, please let me know_
