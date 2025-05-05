%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% local_dirs_mvar_rand_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script must be run in order to run either of those scripts:
% - mec_mvar_rand_cross_coup_rmi.m
% - mec_mvar_rand_cross_coup_rmi_plotting.m
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DIRECTORIES FOR MVAR MODELS
basedir = ['/media/nadinespy/NewVolume1/work/phd/projects/' ...
	'mec_mvar_km/mec_simulations'];

addpath([basedir '/code/multivar_autoreg/rand_cross_coup_rmi/']);
addpath([basedir '/code/functions']);

%% AMEND STRING VARS (FOR COMPUTING AND/OR PLOTTING)

% define array of measure names as used in simulations for plotting 
% & run necessary parameter scripts

% full set:
params_dd_ce_co_info_mvar_rand_cross_coup_rmi
params_integrated_info_mvar_rand_cross_coup_rmi
params_control_mvar_rand_cross_coup_rmi
params_phiid_mmi_mvar_rand_cross_coup_rmi
params_phiid_ccs_mvar_rand_cross_coup_rmi

% measure_in_simulation = {'phiid_mmi', 'phiid_ccs'}; 
% , 'integrated_info', 'dd_ce_co_info' , 'control'};

measure_in_simulation = {'phiid_ccs', 'integrated_info'};

%{
% RUN THIS ONLY WHEN COMPUTING
% measure (not as in the Matlab structs, but as used in the 
% bash/Matlab script for computing)
measure_in_filename = 'dd_ce_co_info';

% define 'measure' for mec_mvar_rand_cross_coup_rmi.m
measure = measure_in_filename;
%}

%% PLOTTING/DATA-SPECIFIC DIRECTORIES

pathout_plots_measures	= [basedir '/results/plots/' ...
	'multivar_autoreg/rand_cross_coup_rmi/' ...
	model_specific_path '/'];
pathout_data_measures	= [basedir '/results/analyses/' ...
	'multivar_autoreg/rand_cross_coup_rmi/' ...
	model_specific_path '/'];

%% FOR PLOTTING

% norm values (exponential scale)
all_cross_coup = 2.^linspace(cross_coup_range(1), ...
	cross_coup_range(2), n_cross_coup); 
% rmi values (exponential scale)
all_rmi   = 2.^linspace(rmi_range(1), rmi_range(2), n_rmi); 

% choose step size for displaying x- and y-tick values 
rmi_ticks_steps = n_rmi/10;
cross_coup_ticks_steps = n_cross_coup/10;


