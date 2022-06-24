%% PLOTTING FOR KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements practical synergy capacity for 256-node Kuramoto oscillators with different couplings and phase lags, 
% two different top-level macro variables (variance of synchronies & global average pairwise synchrony between communities), 
% two mid-level macro variables (synchronies & pairwise synchronies), and four different micro variables (phases, raw signal, 
% synchronies, and pairwise synchronies).

% major sections in this script:
%	- distributions of micro and macro variables						--> get_all_kuramoto_distr_plots()	
%
%	- heatmaps for correlations between micro & macro variables				--> get_all_kuramoto_corr_heatmaps();
%														    uses plot_heatmaps()
%
%	- heatmaps for practical CE, DC, and CD							--> get_all_kuramoto_practCE_heatmaps();
%														    uses plot_heatmaps()
%
%	- heatmaps for dynamical independence							--> get_all_kuramoto_DD_heatmaps();
%														    uses plot_heatmaps()
%
%	- scatter plots for practical CE, with fixed A, and varying beta			--> get_all_kuramoto_practCE_scatterplots();
%														    uses plot_scatterplos_measure()
%
%	- scatter plots for dynamical dependence, with fixed A, and varying beta	--> get_all_kuramoto_DD_scatterplots();
%														    uses plot_scatterplos_measure()
%		
%	- scatter plots for sigma met mean & sigma chi mean,					--> get_all_kuramoto_met_chi_scatterplots();
%       with fixed A, and varying beta								    uses plot_scatterplot_met_chi()

clear all;
clear java;
close all;
clc;

% initialize necessary paths and directories for generated data & plots
n_oscillators = '12';
get_all_kuramoto_directories(); 

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
																				
all_npoints =			[10000]; %, 10000];		% number of data points in time-series
taus =				[3]; %, 10, 100];			% time-lags

%% plotting parameter specification

% data and number of bins for plotting histograms of micro and macro variables
% data can be 'log' or 'raw'
nbins = 100;								
data = 'raw'; 

% axes ticks for heatmaps with beta on x-axis, and A on y-axis

x_axis_heatmaps = {'0.04', '', '', '0.16', '', '', '0.28', '', '', '0.4'};
y_axis_heatmaps = {'0.08', '', '', '0.32', '', '', '0.56', '', '', '0.8'};

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

%% plotting
%% distribution plots of micro and macro variables (selected nodes)

% loop over all values of npoints

% {
  
get_all_kuramoto_distr_plots(data, nbins, network, A_vec, beta_vec, all_npoints, ...
		pathout_data_sim_time_series, pathout_plots_distributions);
		
%}

%% heatmaps for metastability & chimera index

% {

get_all_kuramoto_met_chi_heatmaps(network, all_npoints, x_axis_heatmaps, y_axis_heatmaps, ...
	x_label_heatmaps, y_label_heatmaps, pathout_data_sync, pathout_plots_sigma_chi, ...
	pathout_plots_sigma_met)

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
	x_label_heatmaps, y_label_heatmaps, pathout_data_pract_ce, pathout_plots_pract_ce)

%}

%% dynamical dependence for non-Gaussian continuous variables

% loop over all values of npoints, tau, A, beta

% {

get_all_kuramoto_DD_heatmaps(network, taus, all_npoints, x_axis_heatmaps, y_axis_heatmaps, x_label_heatmaps, ...
	y_label_heatmaps, pathout_data_dd, pathout_plots_dd_mean_pair_sync, pathout_plots_dd_sigma_chi, ...
	pathout_plots_dd_pair_sync, pathout_plots_dd_synchrony, pathout_plots_dd_control_vars)

%}
	
%% scatter plots for practical CE for discretized variables, with fixed A, and varying beta

% loop over all values of npoints, tau

% {

get_all_kuramoto_practCE_scatterplots(network, A_vec, beta_vec, taus, all_npoints, x_label_practCE_scatterplots, ...
		y_label_practCE_scatterplots, pathout_data_pract_ce, pathout_plots_pract_ce_sigma_chi, ...
		pathout_plots_pract_ce_mean_pair_sync, pathout_plots_pract_ce_control_vars)
	
%} 

%% scatter plots for dynamical dependence with non-Gaussian, continuous variables, with fixed A, and varying beta

% loop over all values of npoints, tau

% {

get_all_kuramoto_DD_scatterplots(network, A_vec, beta_vec, taus, all_npoints, x_label_DD_scatterplots, ...
		y_label_DD_scatterplots, pathout_data_dd, pathout_plots_dd_mean_pair_sync, pathout_plots_dd_sigma_chi, ...
		pathout_plots_dd_pair_sync, pathout_plots_dd_synchrony, pathout_plots_pract_ce_control_vars)
	
%} 

%% scatter plots for sigma met mean & sigma chi mean, with fixed A, and varying beta

% loop over all values of npoints 

% {

get_all_kuramoto_met_chi_scatterplots(network, A_vec, beta_vec, all_npoints, x_label_chi_met, y_label_chi_met, ...
		pathout_data_sync, pathout_plots_sigma_chi, pathout_plots_sigma_met)
	
%}
