%% TO DO

% - add integrated information measures?
%
% - generating/defining micro and macro variables only once at the		DONE (but could be better, as I currently pass on
% beginning, and keeping all other scripts flexible for varying			a set of variables to every function, even if some 
% numbers/names of micro and macro variables						variables of that set are not needed in some of them)
%  
% - yeat another step would be to keep model parameters flexible,			DONE (also for measure parameters)
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
% - model parameters (e. g., values for noise correlation & coupling 
% matrix) for CE, DD should be explicit in structs

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

%% create coupling matrices

% {

coupling_matrices = get_coupling_matrices(get_coupling_matrix, coupling_params);
									
%}

%% simulate Kuramoto oscillators, and calculate all micro and macro variables - for all values of model_param1 and model_param2

% GET MICRO AND MACRO VARIABLES (variable names in brackets 
% - to be used to define variable names below):

%	- phases					(phase)		(MICRO)
%	- synchronies				(sync)		(MICRO)
%	- raw signal (cos(phase))		(raw)			(MICRO)
%	- RICA of phases (6 dims)		(rica6_phase)	(MICRO)
%	- RICA of phases (12 dims)		(rica12_phase)	(MICRO)

%	- average pairwise synchrony		(mp_sync)		(MACRO)
%	- chimera-index				(chi)			(MACRO)
%	- sum of RICA of phases (6 dims)	(sum_rica6_phase)	(MACRO)
%	- sum of RICA of phases (12 dims)	(sum_rica12_phase)(MACRO)
%	- sum of phases				(sum_phase)		(MACRO)


% {

% get_all_kuramoto_variables() is the script to generate all micro and 
% macro variables; modify according to what micro and macro variables you 
% wish to have

% get_all_kuramoto_variables() will generate 'phase', 'raw', 'sync', 
% p_sync', 'mp_sync' and 'chi' as variables

get_variables(network, npoints, model_params, coupling_matrices, ...
		pathout_data_sim_time_series, pathout_data_sync);

%}

%% define what micro and macro variables to use to calculate emergence measures

% {

% group names of variables generated in get_all_variables() into micro and macro variabels
micro_variable_names = {'raw', 'phase', 'sync', 'rica6_phase', 'rica12_phase'};
macro_variable_names = {'mp_sync', 'chi', 'sum_phase', 'sum_rica6_phase', 'sum_rica12_phase'};

%}

%% load all micro and macro variables

%{

[phase, chi, sync, p_sync, mp_sync, raw] = load_all_kuramoto_variables(network, model_param1, model_param2, npoints, pathout_data_sim_time_series, ...
	pathout_data_sync);

%}

%% get mean of variance of synchronies across time & mean of variance of synchronies across communities

% {

get_km_met_chi(network, npoints, model_params, ...
		pathout_data_sync)
	
%}

%% quantilize micro and macro variables 

% { 

% generated file names consist of 
% [n_oscillators + network + '_quant' + bin + micro/macro variable name + measure_params{1} + measure_params{2} + time_series_length]
get_all_quant_variables(network, npoints, measure_params, model_params, ...
	micro_variable_names, macro_variable_names, pathout_data_sim_time_series);

% discretizing
% get_all_disc_variables(network, model_param1, model_param2, npoints, micro_variable_names, ...
% 	macro_variable_names, bin_number, pathout_data_sim_time_series);

% binarizing
% get_all_kuramoto_bin_variables(network, model_param1, model_param2, npoints, bin_threshold_phase, ...
% 		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
% 		bin_threshold_sigma_chi, pathout_data_sim_time_series);

%}

%% log transform micro and macro variables

%{

get_all_kuramoto_log_variables(network, model_param1, model_param2, npoints, pathout_data_sim_time_series)

%}

%% calculate average covariance & average correlation between micro & macro variables

% loop over all values of model_param1, model_param2, npoints

%{

get_all_kuramoto_mean_covcorr(network, model_param1, model_param2, npoints, ...
	pathout_data_sim_time_series, ...
	pathout_data_mean_corr, ...
	pathout_data_mean_cov);

%}
	
%% calculate practical CE for quantilized variables

% loop over all values of npoints, measure_param1, model_param1, model_param2

% {

% choose specific file name which will consist of 
% [n_oscillators + network + '_all_practCE_' + ce_veriable_name + time_series_length + tau]

ce_variable_name = 'blubb';

get_all_practCE(network, npoints, measure_params, model_params, micro_variable_names, ...
	macro_variable_names, ce_variable_name, pathout_data_sim_time_series, pathout_data_pract_ce);

%}

%% calculate dynamical dependence for quantilized variables

% {

% choose specific file name which will consist of 
% [n_oscillators + network + '_all_practCE_' + ce_veriable_name + time_series_length + tau]
dd_variable_name = 'blubb';

get_all_DD(network, npoints, measure_params, measure_params_dd, model_params, micro_variable_names, ...
	macro_variable_names, dd_variable_name, pathout_data_sim_time_series, pathout_data_dd);

%}
