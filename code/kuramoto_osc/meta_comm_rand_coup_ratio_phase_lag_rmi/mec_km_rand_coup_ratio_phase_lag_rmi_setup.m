%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mec_km_rand_coup_ratio_phase_lag_rmi_local_dirs.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script must be run in order to run either of those scripts:
% - mec_km_rand_coup_ratio_phase_lag_rmi.m
% - mec_km_rand_coup_ratio_phase_lag_rmi_plotting.m
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DIRECTORIES FOR KM MODELS
basedir = ['/media/nadinespy/NewVolume1/work/phd/projects/' ...
	'mec_mvar_km/mec_simulations'];

addpath([basedir '/code/functions']);
addpath([basedir '/code/kuramoto_osc/' ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi']);

%% AMEND STRING VARS (FOR COMPUTING AND/OR PLOTTING)

% define array of measure names as used in simulations for plotting
% and run parameter script

% full set:

params_dd_ce_co_info_km_rand_coup_ratio_phase_lag_rmi
measure_in_simulation = {'dd_ce_co_info'};

%{
% RUN THIS ONLY WHEN COMPUTING
% measure (not as in the Matlab structs, but as used in the 
% bash/Matlab script for computing)
measure_in_filename = 'dd_ce_co_info';

% define 'measure' for mec_km_rand_coup_ratio_phase_lag_rmi.m
measure = measure_in_filename;
%}

%% PLOTTING/DATA-SPECIFIC DIRECTORIES
pathout_plots_measures	= [basedir '/results/plots/kuramoto_osc/' ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi/' ...
	model_specific_path '/'];
pathout_data_measures	= [basedir '/results/analyses/' ...
	'kuramoto_osc/meta_comm_rand_coup_ratio_phase_lag_rmi/' ...
	model_specific_path '/'];

%% FOR PLOTTING

% rmi values (exponential scale)
all_rmi = 2.^linspace(rmi_range(1), rmi_range(2), n_rmi);
% phase lag - flip order so that values increase
all_phase_lag = flip(pi/2 - linspace(phase_lag_range(1), ...
	phase_lag_range(2), n_phase_lag));
% coupling ratio mean values
all_coup_ratio_mean = linspace(coup_ratio_mean_range(1), ...	
	coup_ratio_mean_range(2), n_coup_ratio_mean);

% choose what number of different coupling ratio means to plot 
% (e. g., if choosing to plot 50 different heatmaps for 50
% different coup ratio mean values, then 50 coupling ratio mean 
% values are evenly sampled from the whole range of the original 
% full set of coupling ratio mean values (which might, say, 100))
n_rmi_to_plot = n_rmi;
sample_rmi = linspace(rmi_range(1), ...
	rmi_range(2), n_rmi_to_plot);

% choose step size for displaying x- and y-tick values 
coup_ratio_mean_ticks_steps = n_coup_ratio_mean/10;
phase_lag_ticks_steps = n_phase_lag/10;
