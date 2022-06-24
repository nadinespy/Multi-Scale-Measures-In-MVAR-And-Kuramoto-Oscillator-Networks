%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?
% - add integrated information measures?

% - generating/defining micro and macro variables only once at the beginning (maybe hard-coding them only in one script), 
%   and keeping all other scripts flexible for varying numbers/names of micro and macro variables would be a next thing to do
% - yeat another step would be to keep model parameters flexible, and only hard-code them at the beginning

% - add Kraskov estimation for practical CE, and see whether multi-dimensional macros are doable
% - search for best subsets in micro and macro that maximize DD and practical CE

% - for Kuramoto oscillators with lots of nodes: do dim reduction
% - store time-series, CE, DD etc. in structs for all model parameters?

%% KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

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

cd /media/nadinespy/NewVolume1/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/EmergenceComplexityMeasuresComparison_Matlab/scripts

% initialize necessary paths and directories for generated data & plots
n_oscillators = '12';
get_all_kuramoto_directories(); 

%% model parameter specification 

% model name for saving files
network = '12kuramoto';

% kuramoto parameter specification
intra_comm_size =			4;					% intra-community size
n_communities =			3;					% number of communities		
A_vec =				linspace(0.08, 0.8, 10);	% A vector
beta_vec =				linspace(0.04, 0.4, 10);	% noise correlation vector: use beta values only up to 0.4, as sigma met & 
										% sigma chi turn out to be zero for greater values of beta; in these cases, 
										% sigma chi will be a non-varying zero macro variable, yielding erroneous 
										% values for emergence 
									
all_npoints =			[10000]; %, 10000];		% number of data points in time-series
taus =				[3]; % , 10, 100];		% time-lags

% parameters for different discretization methods
disc_method =			['quant']; %, 'disc', 'bin'];
bin_number =			[4]; % , 3, 5, 7];		% number of bins to do discretization for method 'disc'
quantile_number =			[5]; % , 3, 4, 7];		% number of quantiles to do discretization for method 'quant'

% thresholds to do discretization for method 'bin'
bin_threshold_phase =		0;
bin_threshold_raw_signal =	0;
bin_threshold_sync =		0.85;
bin_threshold_pair_sync =	0.85;
bin_threshold_sigma_chi =	0.15;

%% measure parameter specification

% method for practical CE - 'Gaussian' or 'discrete'
method_practCE =	['discrete']; % , 'Gaussian']; % to be expanded with 'Kraskov'

% method for dynamical dependence (DD) - 'Kraskov', 'Gaussian', or 'Discrete'
method_DD =		['Kraskov']; % 'Gaussian', 'Discrete'];
tau_steps =		[1]; % , 2, 3];
kraskov_param =	[4]; % , 5, 6];

%% create coupling matrices for each A

% {

% get coupling matrices
for i = 1:length(A_vec)
	A = A_vec(i);	
	kuramoto_coupling_matrices(:,:,i) = get_kuramoto_coupling_matrix(intra_comm_size, n_communities, A);
end 
										
%}

%% simulate Kuramoto oscillators, and calculate all micro and macro variables - for all values of A and all values of beta

% GET MICRO AND MACRO VARIABLES:
%	- phases				(MICRO)
%	- synchronies			(MICRO)
%	- raw signal (cos(phase))	(MICRO)
%	- average pairwise synchrony	(MACRO)
%	- chimera-index			(MACRO)

% {

get_all_kuramoto_variables(network, intra_comm_size, n_communities, kuramoto_coupling_matrices, ...
		all_npoints, A_vec, beta_vec, pathout_data_sim_time_series, pathout_data_sync);

%}

%% load all micro and macro variables

%{

[phase, sigma_chi, synchrony, pair_sync, mean_pair_sync, raw_signal, shuffled_sigma_chi, shuffled_phase] = load_all_kuramoto_variables(network, A_vec, beta_vec, all_npoints, pathout_data_sim_time_series, ...
	pathout_data_sync);

%}

%% get mean of variance of synchronies across time & mean of variance of synchronies across communities

% {

get_all_kuramoto_met_chi(network, A_vec, all_npoints, ...
		pathout_data_sync)
	
%}

%% quantilize micro and macro variables 

% { 

micro_variable_names = {'raw_signal', 'phase', 'synchrony'};
macro_variable_names = {'mean_pair_sync', 'sigma_chi'};

get_all_quant_variables(network, all_npoints, A_vec, beta_vec, micro_variable_names, ...
	macro_variable_names, quantile_number, pathout_data_sim_time_series);

% discretizing
% get_all_disc_variables(network, A_vec, beta_vec, all_npoints, micro_variable_names, ...
% 	macro_variable_names, bin_number, pathout_data_sim_time_series);

% binarizing
% get_all_kuramoto_bin_variables(network, A_vec, beta_vec, all_npoints, bin_threshold_phase, ...
% 		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
% 		bin_threshold_sigma_chi, pathout_data_sim_time_series);

%}

%% log transform micro and macro variables

%{

get_all_kuramoto_log_variables(network, A_vec, beta_vec, all_npoints, pathout_data_sim_time_series)

%}

%% calculate average covariance & average correlation between micro & macro variables

% loop over all values of A, beta, npoints

%{

get_all_kuramoto_mean_covcorr(network, A_vec, beta_vec, all_npoints, ...
	pathout_data_sim_time_series, ...
	pathout_data_mean_corr, ...
	pathout_data_mean_cov);

%}
	
%% calculate practical CE for quantilized variables

% loop over all values of npoints, tau, A, beta

% {

get_all_kuramoto_practCE(network, all_npoints, taus, A_vec, beta_vec, method_practCE, disc_method, ...
		micro_variable_names, macro_variable_names, pathout_data_sim_time_series, pathout_data_pract_ce)
	
%}

%% calculate dynamical dependence for non-Gaussian continuous variables

% {

get_all_kuramoto_DD(network, A_vec, beta_vec, kuramoto_coupling_matrices, taus, tau_steps, all_npoints, ...
		method_DD, pathout_data_sim_time_series, pathout_data_dd, kraskov_param)

%}
