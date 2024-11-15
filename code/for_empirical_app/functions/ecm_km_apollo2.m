%% KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements practical synergy capacity for 256-node Kuramoto oscillators with different couplings and phase lags, two different macro variables 
% (variance of synchronies & global average pairwise synchrony between communities), and three different micro variables (thetas, cos(thetas), synchronies, and 
% binarized synchronies).

% major sections in this script:
%	- choice of parameters (time-lag, data length, and thresholds)
%put 
%	- create coupling matrices & noise correlation vectors					--> get_kuramoto_coupling_matrix()
%
%	- simulate Kuramoto models, including all micro and macro variables,			--> get_all_kuramoto_variables();
%	  for all values of A and all values of beta							    uses get_kuramoto_mean_pair_sync(),
%															    sim_kuramoto_oscillators(), and
%															    shuffle_rows()
%
%	- binarize micro and macro variables								--> get_all_kuramoto_bin_variables();
%															    uses get_binarized_variables()
%
%	- log transform micro and macro variables								--> get_all_kuramoto_log_variables()
%
%	- calculate average covariance & average correlation of					--> get_all_kuramoto_mean_covcorr();
%	  micro & macro variables										    uses get_mean_covcorr_one_matrix() and
%															    get_mean_covcorr_two_matrices()
%
%	- calculate practical CE										--> get_all_kuramoto_practCE();
%															    uses get_practCE()
%
%	- calculate dynamical indepndence									--> get_all_kuramoto_DD();
%															    uses get_DD()
%	- plotting
%		- distributions of micro and macro variables						--> get_all_kuramoto_distr_plots()	
%
%		- heatmaps for correlations between micro & macro variables				--> get_all_kuramoto_corr_heatmaps();
%															    uses plot_heatmaps()
%
%		- heatmaps for practical CE, DC, and CD							--> get_all_kuramoto_practCE_heatmaps();
%															    uses plot_heatmaps()
%
%		- heatmaps for dynamical independence							--> get_all_kuramoto_DD_heatmaps();
%															    uses plot_heatmaps()
%
%		- scatter plots for practical CE, with fixed A, and varying beta			--> get_all_kuramoto_practCE_scatterplots();
%															    uses plot_scatterplos_measure()
%
%		- scatter plots for dynamical dependence, with fixed A, and varying beta	--> get_all_kuramoto_DD_scatterplots();
%															    uses plot_scatterplos_measure()
%		
%		- scatter plots for sigma met mean & sigma chi mean,					--> get_all_kuramoto_met_chi_scatterplots();
%		  with fixed A, and varying beta								    uses plot_scatterplot_met_chi()

clear all;
clear java;
close all;
clc;

cd '/its/home/ns508/job_templates/batch_serial'
addpath '/its/home/ns508/job_templates/batch_serial/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

% directories for generated data
pathout_data = ['/its/home/ns508/job_templates/batch_serial/output/analyses/'];

%% model parameter specification 

% model name for saving files
network =				'256km';

% kuramoto parameter specification
intra_comm_size =			32;					% intra-community size
n_communities =			8;					% number of communities		
A =					linspace(0.08, 0.8, 100);	% vector with different values for A
beta =				linspace(0.04, 0.4, 100);	% vector with different values for noise correlation: use beta values only 
										% up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
										% values of beta; in these cases, sigma chi will be a non-varying 
										% zero macro variable, yielding erroneous values for emergence 
										
npoints =				[2000, 10000, 100000]; %, 10000];	% number of data points in time-series

% parameters for different discretization methods
disc_method =			{'quant'}; %, 'even', 'bin']; % choose discretization method: 'quant' for using quantiles, 
										% 'even' for discretizing into evenly spaced sections of the state space, 
										% 'bin' for binarizing (scripts for latter two not yet modified)
										
bins =				[1, 3, 7]; % , 3, 4, 7];% number of bins to do discretization for method 'quant' and 'disc'

% % thresholds to do discretization for method 'bin'		% not recently modified
% bin_threshold_phase =		0;
% bin_threshold_raw =		0;
% bin_threshold_sync =		0.85;
% bin_threshold_p_sync =	0.85;
% bin_threshold_chi =		0.15;

%% measure parameter specification

tau =					[1, 3, 100]; % , 10, 100];% time-lags

% method for measures based on standard Shannon-information (i.e., Dynamical Independence (DD), and 
% practical Causal Emergence (practCE)); can be 'Discrete', 'Gaussian' or 'Kraskov' 
% ('Kraskov' so far only works for dynamical independence)
method_standard_mi =		{'Discrete'}; % , 'Gaussian', 'Kraskov'}; % to be expanded with 'Kraskov' for practCE
kraskov_param =			[2, 3, 4]; % , 5, 6];	% not yet implemented in loops

% not yet built into the loops
tau_steps =				[1, 3, 10]; % , 2, 3];

%% assign generic variable names further used in the script

% modify the following functions according to model
get_variables =			@get_kuramoto_variables;
get_coupling_matrix =		@get_kuramoto_coupling_matrix;
get_coupling_matrices =		@get_kuramoto_coupling_matrices;

% put all coupling parameters into one cell structure
coupling_params =			{A, intra_comm_size, n_communities};		% must be in that order for 

% put all model parameters into one cell structure
model_params =			{A, beta, intra_comm_size, n_communities};	% model parameters for kuramoto oscillators;
													% must be in that order

% put all measure parameters common to practCE & DD into one cell structure
measure_params =			{tau, method_standard_mi, kraskov_param, ...	% measure parameters common to both practCE and DD;
					disc_method, bins};					% must be in that order

% put measure parameters specific to DD into one cell structure
measure_params_dd =		{tau_steps};						

%% print out current job

fprintf('this job does a parameter sweep over: npoints: %d, tau: %d, method_standard_mi: %d, kraskov_param: %d, disc_method: %d, bins: %d, tau_steps: %d\n', npoints, tau, method_standard_mi, kraskov_param, disc_method, bins, tau_steps);

%% create coupling matrices

% {

coupling_matrices = get_coupling_matrices(get_coupling_matrix, coupling_params);
									
%}


%% simulate Kuramoto oscillators, and calculate all micro and macro variables - for all values of model_param1 and all values of model_param2

% GET MICRO AND MACRO VARIABLES (variable names in brackets):
%	- phases				(phase)	(MICRO)
%	- synchronies			(sync)	(MICRO)
%	- raw signal (cos(phase))	(raw)		(MICRO)
%	- average pairwise synchrony	(mp_sync)	(MACRO)
%	- chimera-index			(chi)		(MACRO)

% {

% get_all_kuramoto_variables() is the script to generate all micro and macro variables; modify according to what
% micro and macro variables you wish to have

% get_all_kuramoto_variables() will generate 'phase', 'raw', 'sync', p_sync', 'mp_sync' and 'chi' as variables
get_variables(network, model_params, coupling_matrices, npoints, ...
	pathout_data, pathout_data);

% group names of variables generated in get_all_variables() into micro and macro variabels
micro_variable_names = {'raw', 'phase', 'sync'};
macro_variable_names = {'mp_sync', 'chi'};

%}

%% quantilize micro and macro variables 

% { 

get_all_quant_variables(network, npoints, measure_params, model_params, ...
	micro_variable_names, macro_variable_names, pathout_data);

% discretizing
% get_all_disc_variables(network, model_param1, model_param2, npoints, micro_variable_names, ...
% 	macro_variable_names, bin_number, pathout_data_sim_time_series);

% binarizing
% get_all_kuramoto_bin_variables(network, model_param1, model_param2, npoints, bin_threshold_phase, ...
% 		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
% 		bin_threshold_sigma_chi, pathout_data_sim_time_series);

%}

%% binarize micro and macro variables 

%{ 

get_all_kuramoto_bin_variables(network, A_vec, beta_vec, all_npoints, bin_threshold_phase, ...
		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
		bin_threshold_sigma_chi, pathout_data);

%}

%% log transform micro and macro variables

%{

get_all_kuramoto_log_variables(network, A_vec, beta_vec, all_npoints, pathout_data)

%}

%% calculate average covariance & average correlation between micro & macro variables

% loop over all values of A, beta, npoints

%{

get_all_kuramoto_mean_covcorr(network, A_vec, beta_vec, all_npoints, ...
	pathout_data, ...
	pathout_data, ...
	pathout_data);

%}
	
%% calculate practical CE for quantilized variables

% loop over all values of npoints, measure_param1, model_param1, model_param2

% {

get_all_practCE(network, npoints, measure_params, model_params, micro_variable_names, ...
	macro_variable_names, pathout_data, pathout_data)
	
%}

%% calculate dynamical dependence for quantilized variables

% {

get_all_DD(network, npoints, measure_params, measure_params_dd, model_params, ...
	micro_variable_names, macro_variable_names, pathout_data, pathout_data)

%}
