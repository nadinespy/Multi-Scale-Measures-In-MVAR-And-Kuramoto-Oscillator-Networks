%% TO DO

% - add integrated information measures?
%
% - generating/defining micro and macro variables only once at the		DONE (but could be better, as I currently pass on
% beginning, and keeping all other scripts flexible for varying			a set of variables to every function, even if some 
% numbers/names of micro and macro variables						variables of that set are not needed in some of them)
%  
% - yet another step would be to keep model parameters flexible,			DONE (also for measure parameters)
% and only hard-code them at the beginning 
%
% - get model parameters for which variables are Gaussian, and get DD/CE
% only for those
%
% - add Kraskov estimation for practical CE, and see whether 
% multi-dimensional macros are doable
%
% - add discrete case for DD using JIDT
%
% - search for best subsets in micro and macro that maximize DD and CE
%
% - for Kuramoto oscillators with lots of nodes: do dim reduction
%
% - model parameters (e. g., values for noise correlation & coupling 
% matrix) for CE, DD should be explicit in structs

%% KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% DESSCRIPTION NEEDS TO BE UPDATED

% This script implements practical synergy capacity for 256-node Kuramoto oscillators with different couplings and phase lags, 
% two different top-level macro variables (variance of synchronies & global average pairwise synchrony between communities), 
% two mid-level macro variables (synchronies & pairwise synchronies), and four different micro variables (phases, raw signal, 
% synchronies, and pairwise synchronies).

% major sections in this script:
%	- choice of parameters (time-lag, data length, and thresholds)
%
%	- create coupling matrices & noise correlation vectors					--> get_kuramoto_coupling_matrix()
%
%	- simulate Kuramoto models, including all micro and macro variables,			--> get_all_kuramoto_variables();
%	  for all values of model_param1 and all values of model_param2							    uses get_kuramoto_mean_pair_sync(),
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
%		- scatter plots for practical CE, with fixed model_param1, and varying model_param2			--> get_all_kuramoto_practCE_scatterplots();
%															    uses plot_scatterplos_measure()
%
%		- scatter plots for dynamical dependence, with fixed model_param1, and varying model_param2	--> get_all_kuramoto_DD_scatterplots();
%															    uses plot_scatterplos_measure()
%		
%		- scatter plots for sigma met mean & sigma chi mean,					--> get_all_kuramoto_met_chi_scatterplots();
%		  with fixed model_param1, and varying model_param2								    uses plot_scatterplot_met_chi()

clear all;
clear java;
close all;
clc;

cd /media/nadinespy/NewVolume1/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/EmergenceComplexityMeasuresComparison_Matlab/scripts
addpath(genpath('/media/nadinespy/NewVolume1/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/kuramoto'))

directories = @get_all_kuramoto_directories;

% initialize necessary paths and directories for generated data & plots
n_oscillators = '12';
directories(); 

%% model parameter specification 

% model name for saving files (type and size of model)
network =				'12km';

% kuramoto parameter specification
intra_comm_size =			4;					% intra-community size
n_communities =			3;					% number of communities		
A =					linspace(0.08, 0.8, 10);		% vector with different values for A
beta =				linspace(0.08, 0.8, 10);		% vector with different values for noise correlation: use beta values only 
										% up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
										% values of beta; in these cases, sigma chi will be a non-varying 
										% zero macro variable, yielding erroneous values for emergence 
										
npoints =				[10000]; %, 10000];		% number of data points in time-series

%% measure parameter specification

% parameters for different discretization methods
disc_method =			{'quant'}; %, 'even', 'bin']; % choose discretization method: 'quant' for using quantiles, 
										% 'even' for discretizing into evenly spaced sections of the state space, 
										% 'bin' for binarizing (scripts for latter two not yet modified)
										
bins =				[1, 3, 7]; % , 1, 7];		% number of bins to do discretization for method 'quant' and 'disc'

% thresholds to do discretization for method 'bin'		% not recently modified
% bin_threshold_phase =		0;
% bin_threshold_raw =		0;
% bin_threshold_sync =		0.85;
% bin_threshold_p_sync =	0.85;
% bin_threshold_chi =		0.15;


tau =					[1, 3, 10]; % , 1, 10];		% time-lags

% method for measures based on standard Shannon-information (i.e., Dynamical
% Independence (DD), and practical Causal Emergence (practCE)); can be 
% 'Discrete', 'Gaussian' or 'Kraskov' ('Kraskov' so far only works for DD)

method_standard_mi =		{'Gaussian'};		% , 'Gaussian', 'Kraskov', 'Discrete'}; % to be expanded with 'Kraskov' for practCE
method_pid_mi =			{'MMI', 'CCS'};
kraskov_param =			[2, 3, 4];			% , 3, 4];	% not yet implemented in loops

% not yet built into the loops
tau_steps =				[1, 3, 10]; % , 3, 10];

%% assign generic variable names further used in the script

% modify the following functions according to model
get_variables =			@get_km_variables;					% function to get micro/macro variables of interest
get_coupling_matrix =		@get_km_coupling_matrix;				% function to get coupling matrix for one model
get_coupling_matrices =		@get_km_coupling_matrices;				% function to get coupling matrices for all models

% put all coupling parameters into one cell structure
coupling_params =			{A, intra_comm_size, n_communities};		% must be in that order 

% put all model parameters into one cell structure
model_params =			{A, beta, intra_comm_size, n_communities};	% model parameters for kuramoto oscillators;
													% must be in that order

% put all measure parameters common to practCE & DD into one cell structure
measure_params =			{tau, method_standard_mi, kraskov_param, ...	% measure parameters common to both practCE and DD;
					disc_method, bins};					% must be in that order

% put measure parameters specific to DD into one cell structure
measure_params_dd =		{tau_steps};						


%% simulate, calculate, plot

% {

% group names of variables generated in get_all_variables() into micro and macro variabels
micro_variable_names = {'raw', 'phase', 'sync', 'rica6_phase', 'rica12_phase'};
macro_variable_names = {'mp_sync', 'chi', 'sum_phase', 'sum_rica6_phase', 'sum_rica12_phase'};

% file prefixes to distinguish different DD & CE struct files, and not overwrite them
ce_variable_name = 'standard';
dd_variable_name = 'standard';

ecm_kuramoto_sim_calc();

ecm_kuramoto_plot();

%}