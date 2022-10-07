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

% DESCRIPTION NEEDS TO BE UPDATED

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
%	  for all values of model_param1 and all values of model_param2				    uses get_kuramoto_mean_pair_sync(),
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
addpath(genpath('/media/nadinespy/NewVolume1/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/mutual_info_kNN'))

directories = @get_all_kuramoto_directories;

% initialize necessary paths and directories
n_oscillators = '12';
directories(); 

%% model parameter specification 

% model name for saving files (type and size of model)
network			= '12km';

% kuramoto parameter specification
intra_comm_size		= 4;						% intra-community size
n_communities		= 3;						% number of communities		
A				= linspace(0.08, 0.8, 1);		% vector with different values for A
beta				= linspace(0.08, 0.8, 1);		% vector with different values for noise correlation: use beta values only 
										% up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
										% values of beta; in these cases, sigma chi will be a non-varying 
										% zero macro variable, yielding erroneous values for emergence 

disc_methods		= {'quant'}; %, 'even', 'bin'};	% choose discretization method: 'quant' for using quantiles, 
										% 'even' for discretizing into evenly spaced sections of the state space, 
										% 'bin' for binarizing (scripts for latter two not yet modified)
										
bins				= [1]; %, 3, 7];				% number of bins to do discretization for method 'quant' and 'disc'

get_coupling_matrix	= @get_km_coupling_matrix;		% specify function to generate one coupling matrix
get_coupling_matrices	= @get_km_coupling_matrices;		% specify function to generate coupling matrices
get_variables		= @get_km_variables;			% specify function to generate micro and macro variables


%% measure parameter specification
% -------------------------------------------------------------------------
% necessary input arguments

measures			= {'phiidCE', 'phiidDC', 'phiidCD', 'shannonCE', 'shannonDC', 'shannonCD', 'DD'};	%, 'ShannonCE', 'ShannonDC', 'ShannonCD', % emergence measures
						% 'phiidCE', 'phiidDC', 'phiidCD'};					
methods			= {'Kraskov', 'Gaussian', 'Discrete'};				% to be expanded with 'Kraskov' for practCE; discrete method can practically 
													% handle only systems with 10-12 binary variables; if variables are not binary,  
													% then even less (as the joint state-space grows)
time_lags			= [3]; %, 3, 10];							% time-lags
time_lengths		= [10000]; %, 2000];

% -------------------------------------------------------------------------
% optional input arguments (depending on method)

kraskov_params		= [2]; %, 3, 4];							

% -------------------------------------------------------------------------
% input arguments specific to measures

% PhiID-CE
red_funcs			= {'MMI'}; %, 'CCS'};

% DD
time_steps			= [1, 3]; %, 3, 10]; 

%% put all parameters into cell structures

% -------------------------------------------------------------------------
% model parameters

% model parameters for simulating kuramoto oscillators; must be in that order
model_sim_params.A			= A ;		 
model_sim_params.beta			= beta;				
model_sim_params.intra_comm_size	= intra_comm_size;
model_sim_params.n_communities	= n_communities;

% model parameters to calculate emergence for; must be in that order
model_calc_params.A			= A ;		 
model_calc_params.beta			= beta;																					

% -------------------------------------------------------------------------														
% put all common measure parameters common to Shannon CE, PhiID-CE & DD 
% into one cell structure

measure_params.measures			= measures;
measure_params.methods			= methods;
measure_params.time_lags		= time_lags;
measure_params.time_lengths		= time_lengths;
measure_params.kraskov_params		= kraskov_params;
measure_params.disc_methods		= disc_methods;
measure_params.bins			= bins;

% -------------------------------------------------------------------------
% put measure parameters specific to DD into one cell structure
measure_params_dd.time_steps		= time_steps;

% put measure parameters specific to PhiID-CE into one cell structure
measure_params_phiid_ce.red_funcs	= red_funcs;

% -------------------------------------------------------------------------
% pathouts for output
pathout.pathout_data_shannon_ce	= pathout_data_shannon_ce;
pathout.pathout_data_phiid_ce		= pathout_data_phiid_ce;
pathout.pathout_data_dd			= pathout_data_dd;

% pathin
pathin.pathout_data_sim_time_series = pathout_data_sim_time_series;
pathin.pathout_data_sync		= pathout_data_sync;

% -------------------------------------------------------------------------
% group names of variables generated in get_all_variables() into 
% micro and macro variabels

% micro_variable_names = {'raw', 'phase', 'sync', 'rica6_phase', ...
% 	'rica12_phase'};
% macro_variable_names = {'mp_sync', 'chi', 'sum_phase', ...
% 	'sum_rica6_phase', 'sum_rica12_phase'};

micro_variable_names = {'phase'};
macro_variable_names = {'mp_sync', 'chi'};

%% get variables (simulate data, and discretize them)

%{
ecm_get_variables();
%}

%% calculate emergence

% file prefixes to distinguish different struct files, and not overwrite them
variable_name = 'standard';

emergence_results = get_all_emergence(network, model_calc_params, measure_params, ...
		micro_variable_names, macro_variable_names, variable_name, ...
		pathin, pathout, 'measure_params_phiid_ce', measure_params_phiid_ce, ...
		'measure_params_dd', measure_params_dd);

%}