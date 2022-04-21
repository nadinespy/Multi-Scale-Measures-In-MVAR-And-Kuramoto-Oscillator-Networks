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

%% KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements practical synergy capacity for 256-node Kuramoto oscillators with different couplings and phase lags, two different macro variables 
% (variance of synchronies & global average pairwise synchrony between communities), and three different micro variables (thetas, cos(thetas), synchronies, and 
% binarized synchronies).

% major sections in this script:
%	- choice of parameters (time-lag, data length, and thresholds)
%
%	- create coupling matrices & noise correlation vectors					--> get_kuramoto_coupling_matrix()
%
%	- simulate Kuramoto models, including all micro and macro variables,			--> get_all_kuramoto_variables();
%	  for all values of A and all values of beta							    uses get_kuramoto_mean_pair_sync(), and
%															    sim_kuramoto_oscillators()
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

% initialize necessary paths and directories for generated data & plots
get_256kuramoto_directories(); 

%% model parameter specification 

% model name for saving files
network = '12kuramoto';

% parameter specification
intra_comm_size =			4;					% intra-community size
n_communities =			3;					% number of communities		
A_vec =				linspace(0.08, 0.8, 10);	% A vector
beta_vec =				linspace(0.04, 0.4, 10);	% noise correlation vector: use beta values only up to 0.4, as sigma met & 
										% sigma chi turn out to be zero for greater values of beta; in these cases, 
										% sigma chi will be a non-varying zero macro variable, yielding erroneous 
										% values for emergence 
																				
all_npoints =			[2000]; %, 10000];		% number of data points in time-series
taus =				[10]; %, 10, 100];		% time-lags

bin_threshold_phase =		0;
bin_threshold_raw_signal =	0;
bin_threshold_sync =		0.9;
bin_threshold_pair_sync =	0.9;
bin_threshold_sigma_chi =	0.25;

%% measure parameter specification

% method for practical CE - 'Gaussian' or 'discrete'
method_practCE = 'discrete';

% method for dynamical dependence (DD) - 'Kraskov', 'Gaussian', or 'Discrete'
method_DD = 'Kraskov';
tau_steps = 1;
kraskov_param = 4;

%% plotting parameter specification

% data and number of bins for plotting histograms of micro and macro variables
% data can be 'log' or 'raw'
nbins = 100;								
data = 'raw'; 

% axes ticks for heatmaps with beta on x-axis, and A on y-axis

x_axis_heatmaps = {'0.04', '', '', '0.16', '', '', '0.28', '', '', '0.4'};
y_axis_heatmaps = {'0.08', '', '', '0.32', '', '', '0.56', '', '', '0.8'};
%x_axis_heatmaps = {'0.04', '0.4'};
%y_axis_heatmaps = {'0.08', '0.8'};

y_label_heatmaps = 'A';
x_label_heatmaps = 'beta';

% axes ticks for scatterplots with beta on x-axis, and practical CE on y-axis
y_label_practCE_scatterplots = 'practical CE';
x_label_practCE_scatterplots = 'beta';

% axes ticks for scatterplots with beta on x-axis, and dynamical dependence on y-axis
y_label_DD_scatterplots = 'dynamical dependence';
x_label_DD_scatterplots = 'beta';

% axes ticks for scatterplots with beta on x-axis, and sigma chi/sigma met on y-axis
y_label_chi_met = ' ';
x_label_chi_met = 'beta';

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
		A_vec, beta_vec, all_npoints, pathout_data_sim_time_series, pathout_data_sync);

%}

%% binarize micro and macro variables 

% { 

get_all_kuramoto_bin_variables(network, A_vec, beta_vec, all_npoints, bin_threshold_phase, ...
		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
		bin_threshold_sigma_chi, pathout_data_sim_time_series);

%}

%% log transform micro and macro variables

% {

get_all_kuramoto_log_variables(network, A_vec, beta_vec, all_npoints, pathout_data_sim_time_series)

%}

%% calculate average covariance & average correlation between micro & macro variables

% loop over all values of A, beta, npoints

% {

get_all_kuramoto_mean_covcorr(network, A_vec, beta_vec, all_npoints, ...
	pathout_data_sim_time_series, ...
	pathout_data_mean_corr, ...
	pathout_data_mean_cov);

%}
	
%% calculate practical CE for discretized variables

% loop over all values of npoints, tau, A, beta

% {

get_all_kuramoto_practCE(network, A_vec, beta_vec, kuramoto_coupling_matrices, taus, all_npoints, method_practCE, ...
		pathout_data_sim_time_series, pathout_data_pract_ce)

%}

%% calculate dynamical dependence for non-Gaussian continuous variables

% {

get_all_kuramoto_DD(network, A_vec, beta_vec, kuramoto_coupling_matrices, taus, tau_steps, all_npoints, ...
		method_DD, pathout_data_sim_time_series, pathout_data_dd, kraskov_param)

%}

%% plotting

%% distribution plots of micro and macro variables (selected nodes)

% loop over all values of npoints

% {
  
get_all_kuramoto_distr_plots(data, nbins, network, A_vec, beta_vec, all_npoints, ...
		pathout_data_sim_time_series, pathout_plots_distributions);
		
%}

%% correlations between micro & macro variables

% loop over all values of npoints

% {

get_all_kuramoto_corr_heatmaps(network, all_npoints, x_axis_heatmaps, y_axis_heatmaps, ...
	x_label_heatmaps, y_label_heatmaps, pathout_data_mean_corr, pathout_plots_mean_corr)

%}
	
%% practical CE, DC, and CD for discretized variables

% loop over all values of npoints, tau, A, beta

% {

get_all_kuramoto_practCE_heatmaps(network, taus, all_npoints, x_axis_heatmaps, y_axis_heatmaps, ...
	x_label_heatmaps, y_label_heatmaps, pathout_data_pract_ce, pathout_plots_pract_ce_mean_pair_sync, ...
	pathout_plots_pract_ce_sigma_chi)

%}

%% dynamical dependence for non-Gaussian continuous variables

% loop over all values of npoints, tau, A, beta

% {

get_all_kuramoto_DD_heatmaps(network, taus, all_npoints, x_axis_heatmaps, y_axis_heatmaps, x_label_heatmaps, ...
	y_label_heatmaps, pathout_data_dd, pathout_plots_dd_mean_pair_sync, pathout_plots_dd_sigma_chi, ...
	pathout_plots_dd_pair_sync, pathout_plots_dd_synchrony)

%}
	
%% scatter plots for practical CE for discretized variables, with fixed A, and varying beta

% loop over all values of npoints, tau

% {

get_all_kuramoto_practCE_scatterplots(network, A_vec, beta_vec, taus, all_npoints, x_label_practCE_scatterplots, ...
		y_label_practCE_scatterplots, pathout_data_pract_ce, pathout_plots_pract_ce_sigma_chi, ...
		pathout_plots_pract_ce_mean_pair_sync)
	
%} 

%% scatter plots for dynamical dependence with non-Gaussian, continuous variables, with fixed A, and varying beta

% loop over all values of npoints, tau

% {

get_all_kuramoto_DD_scatterplots(network, A_vec, beta_vec, taus, all_npoints, x_label_DD_scatterplots, ...
		y_label_DD_scatterplots, pathout_data_dd, pathout_plots_dd_mean_pair_sync, pathout_plots_dd_sigma_chi, ...
		pathout_plots_dd_pair_sync, pathout_plots_dd_synchrony)
	
%} 

%% scatter plots for sigma met mean & sigma chi mean, with fixed A, and varying beta

% loop over all values of npoints 

% {

get_all_kuramoto_met_chi_scatterplots(network, A_vec, beta_vec, all_npoints, x_label_chi_met, y_label_chi_met, ...
		pathout_data_sync, pathout_plots_sigma_chi, pathout_plots_sigma_met)
	
%}
