%% TO DO

% - put different models into separate scripts
% - add integrated information measures?
% - Shannon-CE for multi-dimensional macros - should already work?
% - integrate PhiIDAnalytical() in get_phiid_atoms()
% - work on the caveats: 
%	- discrete method:
%		- works only for dimensionality-reduced systems when using 
%		  PhiID-CE and DD
%		- can practically handle only systems with 
%		  10-12 binary variables; 
%		- if variables are not binary, then even less 
%		  (as the joint state-space grows)
%	- Kraskov method in PhiID- and Shannon-CE/DC/CD:
%		- works only for systems with even-sized number of variables 
%		  (don’t know why that is the case!)
%		- doesn't work for kraskov in continuous PhiID, as multivariate 
%		  entropy is for normal Gaussians
%     - Shannon-CE/DC/CD: for now, calculations can only consider 1D macro
%       variables

%% MAIN SCRIPT TO SPECIFY CALCULATIONS (LIKE A CONFIG FILE ^^)

% WRITE DESCRIPTION IN MORE DETAIL

% This script implements a specified set of measures of emergence, those measures can be:
%	- Shannon- and PhiID-based Causal Emergence (CE),
%	  Downward Causation (DC), Causal Decoupling (CD),
%	- Dynamical Dependence (DD).
%
% The specified set of measures is calculated for:
%	- continuous data (which can be discretized, too),
%	- a specified range of two model parameters, 
%	- a specified range of common measure parameters,
%	- a specified range of measure parameters that are specific to measures
% 	- specified micro and macro variable combinations,
%	- specified estimation methods (possible ones are Gaussian, Kraskov, Discrete).
%		- Discrete in DD & PhiID-based CE: loads continuous data, does dimensionality reduction 
%		  (two dimensions for micro and one for macro variables), and quantilizes reduced data
%		- Discrete in Shannon-based CE: loads data that have been quantilized before

clear all;
clear java;
close all;
clc;

cd '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/empirical/scripts'
javaaddpath('/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/empirical/functions/infodynamics.jar');
javaaddpath('/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/empirical/scripts/infodynamics.jar');


% get model-specific directories, and determine model size

%**************************************************************************
% KURAMOTO OSCILLATORS
%**************************************************************************
% {
directories		= @get_km_directories;
n_oscillators	= '12';
intra_comm_size	= 4;	% intra-community size
n_communities	= 3;	% number of communities	
%}

%**************************************************************************
% MULTIVARIATE AUTOREGRESSIVE NETWORK
%**************************************************************************
%{
directories = @get_mvar_directories;
n_nodes = '2';
%}

% initialize necessary paths and directories
directories(); 

% -------------------------------------------------------------------------
%% model parameter specification 
% -------------------------------------------------------------------------

%**************************************************************************
% KURAMOTO OSCILLATORS
%**************************************************************************
% {
network			= [n_oscillators 'km'];
	
A				= linspace(0.08, 0.8, 10);		% vector with different values for A
beta				= linspace(0.08, 0.8, 10);		% vector with different values for phase lag 
										% (range 0-0.4: metastable regime; range 0.4-0.8: non-metastable regime)

get_coupling_matrix	= @get_km_coupling_matrix;		% specify function to generate one coupling matrix
get_coupling_matrices	= @get_km_coupling_matrices;		% specify function to generate multiple coupling matrices
get_variables		= @get_km_variables;			% specify function to generate micro and macro variables
%}

%**************************************************************************
% MULTIVARIATE AUTOREGRESSIVE NETWORK
%**************************************************************************
%{
time_lag_for_model	= [1];					% time-lag to simulate MVAR process
network			= [n_nodes 'mvar' ...			% model name
				  '_lag' num2str(time_lag_for_model)];

couplings			= linspace(0.0045, 0.45, 100);	% vector with different coupling values 
noise_corrs			= linspace(0.09, 0.9, 100);		% vector with different values for noise correlation 

get_coupling_matrix	= @get_mvar_coupling_matrix;		% specify function to generate one coupling matrix
get_coupling_matrices	= @get_mvar_coupling_matrices;	% specify function to generate multiple coupling matrices
get_variables		= @get_mvar_variables;			% specify function to generate micro and macro variables
%}

%**************************************************************************
% PARAMETERS COMMON TO ALL MODELS
%**************************************************************************

disc_methods		= {'quant'}; %, 'even', 'bin'};	% choose discretization method: 'quant' for using quantiles, 
										% 'even' for discretizing into evenly spaced sections of the 
										% state space, 'bin' for binarizing (scripts for latter two 
										% not yet modified)
										
bins				= [1]; %, 3, 7];				% number of quantiles or bins to do discretization for methods 
										% 'quant' and 'disc', e. g., using only one quantile would 
										% result in two groups - one below and one above the quantile

% -------------------------------------------------------------------------
%% measure parameter specification
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% necessary input arguments
% -------------------------------------------------------------------------

% caveats: 
% - discrete method:
%		- PhiID-based measures and DD: works only for dimensionality-
%		  reduced systems
%			- can practically handle only systems with 
%			  10-12 binary variables; 
%			- if variables are not binary, then even less 
%		        (as the joint state-space grows)
% - Kraskov method (only applicable to DD & Shannon-CE/DC/CD):
%		- Shannon-CE/DC/CD: works only for systems with even-sized  
%		  number of variables  (don’t know why that is the case!)
%		- PhiID-CE/DC/CD: not implemented yet, as multivariate entropy
%		  is calculated only for Gaussian systems

%measures			= {'DD', 'shannonCE', 'shannonDC', 'shannonCD', ...
%				   'phiidCE', 'phiidDC', 'phiidCD'}; 
measures			= {'phiidCE', 'DD', 'shannonCE'}; 
methods			= {'Discrete', 'Gaussian', 'Kraskov'};%'Kraskov', 'Discrete', 'Gaussian'};
%methods			= {'Discrete'};
time_lags_for_measure	= [1]; %, 3, 10];						
time_lengths		= [2000]; %, 2000];  % for Kuramoto oscillators, this 
							   % results in time_lengths(i)/0.05 
							   % data points for the ith time-length

% -------------------------------------------------------------------------
% optional input arguments (depending on method)
% -------------------------------------------------------------------------

kraskov_params		= [3]; %, 2, 4];							

% -------------------------------------------------------------------------
% input arguments specific to measures
% -------------------------------------------------------------------------

% PhiID-CE
red_funcs			= {'MMI', 'CCS'};%, 'MMI'}; % to be expanded with CCS

% DD
time_steps			= [1]; %, 3, 10]; 

% -------------------------------------------------------------------------
%% put all parameters into cell structures
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% model parameters
% -------------------------------------------------------------------------

%**************************************************************************
% KURAMOTO OSCILLATORS
%**************************************************************************
% {
% model parameters for simulating kuramoto oscillators; must be in that order
model_sim_params.A			= A ;		 
model_sim_params.beta			= beta;				
model_sim_params.intra_comm_size	= intra_comm_size;
model_sim_params.n_communities	= n_communities;

% model parameters to calculate emergence for; must be in that order
model_calc_params.A			= A ;		 
model_calc_params.beta			= beta;	

% pathouts for output
pathout.pathout_data_measures		= pathout_data_measures;
pathout.pathout_data_metastability  = pathout_data_metastability;

% pathin
pathin.pathout_data_sim_time_series = pathout_data_sim_time_series;
pathin.pathout_data_sync		= pathout_data_sync;

% group names of variables generated in get_all_variables() into 
% micro and macro variabels

% micro_variable_names			= {'raw', 'phase', 'community_sync', 'p_sync', ...
% 						      'rica6_phase', 'rica12_phase'};
% macro_variable_names			= {'mp_sync', 'chi', 'sum_phase', 'community_sync', ...
% 							'full_system_sync' 'sum_rica6_phase', ...
% 							'sum_rica12_phase'};

micro_variable_names			= {'phase', 'community_sync'};
macro_variable_names			= {'sum_phase', 'chi', 'mp_sync', 'community_sync', 'full_system_sync'};
																				
%}

%**************************************************************************
% MULTIVARIATE AUTOREGRESSIVE NETWORK 
%**************************************************************************
%{
% model parameters for simulating MVAR networks; must be in that order
model_sim_params.couplings			= couplings;		 
model_sim_params.noise_corrs			= noise_corrs;	
model_sim_params.n_nodes			= str2double(n_nodes);
model_sim_params.time_lag_for_model		= time_lag_for_model;

% model parameters to calculate emergence for; must be in that order
model_calc_params.couplings			= couplings;		 
model_calc_params.noise_corrs			= noise_corrs;	

% pathouts for output
pathout.pathout_data_measures			= pathout_data_measures;

% pathin
pathin.pathout_data_sim_time_series		= pathout_data_sim_time_series;

% group names of variables generated in get_all_variables() into 
% micro and macro variabels
micro_variable_names				= {'nodes'};
macro_variable_names				= {'sum_nodes', 'sum_exp_nodes'};
%}

% -------------------------------------------------------------------------														
% put all measure parameters common to Shannon CE, PhiID-CE & DD 
% into one cell structure
% -------------------------------------------------------------------------

measure_params.measures			  = measures;
measure_params.methods			  = methods;
measure_params.time_lags_for_measure  = time_lags_for_measure;
measure_params.time_lengths		  = time_lengths;
measure_params.kraskov_params		  = kraskov_params;
measure_params.disc_methods		  = disc_methods;
measure_params.bins			  = bins;

% -------------------------------------------------------------------------
% put measure parameters specific to DD into one cell structure
% -------------------------------------------------------------------------

measure_params_dd.time_steps		  = time_steps;

% put measure parameters specific to PhiID-CE into one cell structure
measure_params_phiid_ce.red_funcs	  = red_funcs;

% -------------------------------------------------------------------------
%% get variables (simulate data, and discretize them)
% -------------------------------------------------------------------------

% {

coupling_matrices = get_coupling_matrices(get_coupling_matrix, ...
	model_sim_params);

get_variables(network, ...
	model_sim_params, ...
	measure_params, ...
	coupling_matrices, ...
	pathin);
	
get_all_quant_variables(network, ...
	measure_params, ...
	model_sim_params, ...
	micro_variable_names, ...
	macro_variable_names, ...
	pathin);
	
%}

% -------------------------------------------------------------------------
%% calculate emergence
% -------------------------------------------------------------------------

% file prefixes to distinguish and not overwrite different struct files
struct_prefix = 'all_measures_';

emergence_results = get_all_emergence(network, ...
	model_calc_params, ...
	measure_params, ... 
	micro_variable_names, ...
	macro_variable_names, ...
	pathin, ...
	'measure_params_phiid_ce', measure_params_phiid_ce, ...
	'measure_params_dd', measure_params_dd);

% save emergence_results
save([pathout_data_measures struct_prefix network '_emergence_results.mat'], 'emergence_results'); 

%}


get_km_met_chi_sync(network, model_sim_params, time_lengths, pathin, pathout)