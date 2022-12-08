%% TO DO

% - add integrated information measures?
%
% - get model parameters for which variables are Gaussian, and get DD/CE
% only for those
%
% - Shannon-CE for multi-dimensional macros
%
% - search for best subsets in micro and macro that maximize DD and CE
%
% - for Kuramoto oscillators with lots of nodes: do dim reduction

%% KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% WRITE DESCRIPTION IN MORE DETAIL

% This script implements a specified set of measures of emergence, those measures can be:
%	- Shannon- and PhiID-based Causal Emergence (CE),
%	  Downward Causation (DC), Causal Decoupling (CD),
%	- Dynamical Independence.
%
% The specified set of measures is calculated for:
%	- a specified range of two model parameters, 
%	- a specified range of common measure parameters,
%	- a specified range of measure parameters that are specific to measures
% 	- specified micro and macro variable combinations,
%	- specified estimation methods (possible ones are Gaussian, Kraskov, Discrete).

clear all;
clear java;
close all;
clc;

cd /media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/EmergenceComplexityMeasuresComparison_Matlab/scripts
addpath(genpath('/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/kuramoto'))
addpath(genpath('/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/mutual_info_kNN'))

% kuramoto
%{
directories = @get_km_directories;
n_oscillators = '12';
%}

% mvar
% {
directories = @get_mvar_directories;
n_nodes = '2';
%}

% initialize necessary paths and directories
directories(); 

%% model parameter specification 

% kuramoto
%{
% model name for saving files (type and size of model)
network			= '12km';

intra_comm_size		= 4;						% intra-community size
n_communities		= 3;						% number of communities		
A				= linspace(0.08, 0.8, 10);		% vector with different values for A
beta			= linspace(0.08, 0.8, 10);		% vector with different values for phase lag

get_coupling_matrix	= @get_km_coupling_matrix;		% specify function to generate one coupling matrix
get_coupling_matrices	= @get_km_coupling_matrices;		% specify function to generate multiple coupling matrices
get_variables		= @get_km_variables;			% specify function to generate micro and macro variables
%}

% mvar
% {
% model name for saving files (type and size of model)
network			= '2mvar';

coupling			= linspace(0.08, 0.8, 10);		% vector with different coupling values 
corr_err			= linspace(0.08, 0.8, 10);		% vector with different values for noise correlation 

get_coupling_matrix	= @get_mvar_coupling_matrix;		% specify function to generate one coupling matrix
get_coupling_matrices	= @get_mvar_coupling_matrices;	% specify function to generate multiple coupling matrices
get_variables		= @get_mvar_variables;			% specify function to generate micro and macro variables
%}

disc_methods		= {'quant'}; %, 'even', 'bin'};	% choose discretization method: 'quant' for using quantiles, 
										% 'even' for discretizing into evenly spaced sections of the state space, 
										% 'bin' for binarizing (scripts for latter two not yet modified)
										
bins				= [1, 3]; %, 3, 7];			% number of quantiles or bins to do discretization for methods 'quant' and 'disc', 
										% e. g., using only one quantile would result in two groups - 
										% one below and one above the quantile

%% measure parameter specification
% -------------------------------------------------------------------------
% necessary input arguments

% caveats: 
% - discrete method:
%		- works only for dimensionality reduced systems when using PhiID-CE and DD
%		- can practically handle only systems with 10-12 binary variables; 
%		- if variables are not binary, then even less (as the joint state-space grows)
% - Kraskov method works only for systems with even-sized number of variables
% (don’t know why that is the case!)

measures			= {'shannonCE', 'shannonDC', 'shannonCD','phiidCE', 'phiidDC', 'phiidCD', 'DD'}; 
methods			= {'Kraskov', 'Discrete', 'Gaussian'};
time_lags			= [3, 10]; %, 3, 10];							% time-lags
time_lengths		= [10000]; %, 2000];

% -------------------------------------------------------------------------
% optional input arguments (depending on method)

kraskov_params		= [2, 3]; %, 3, 4];							

% -------------------------------------------------------------------------
% input arguments specific to measures

% PhiID-CE
red_funcs			= {'MMI'}; %, 'CCS'}; % to be expanded with CCS

% DD
time_steps			= [1, 3]; %, 3, 10]; 

%% put all parameters into cell structures

% -------------------------------------------------------------------------
% model parameters

% kuramoto
%{
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

% pathin
pathin.pathout_data_sim_time_series = pathout_data_sim_time_series;
pathin.pathout_data_sync		= pathout_data_sync;

% group names of variables generated in get_all_variables() into 
% micro and macro variabels

% micro_variable_names = {'raw', 'phase', 'sync', 'rica6_phase', ...
% 	'rica12_phase'};
% macro_variable_names = {'mp_sync', 'chi', 'sum_phase', ...
% 	'sum_rica6_phase', 'sum_rica12_phase'};
																				
%}
% -------------------------------------------------------------------------
% mvar 
% {
% model parameters for simulating MVAR networks; must be in that order
model_sim_params.coupling		= coupling ;		 
model_sim_params.corr_err		= corr_err;				

% model parameters to calculate emergence for; must be in that order
model_calc_params.coupling		= coupling ;		 
model_calc_params.corr_err		= corr_err;	

% pathouts for output
pathout.pathout_data_measures		= pathout_data_measures;

% pathin
pathin.pathout_data_sim_time_series = pathout_data_sim_time_series;

% group names of variables generated in get_all_variables() into 
% micro and macro variabels
micro_variable_names = {'nodes' 'rica2_nodes'};
macro_variable_names = {'sum_nodes', 'sum_exp_nodes'};
%}

% -------------------------------------------------------------------------														
% put all measure parameters common to Shannon CE, PhiID-CE & DD 
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

%% get variables (simulate data, and discretize them)

% {

% TO-DO: WRITE COUPLING FUNCTION FOR MVAR
get_variables(network, model_sim_params, measure_params, coupling_matrices, ...
		pathout);
%}

%% calculate emergence

% file prefixes to distinguish and not overwrite different struct files
struct_prefix = '01_';

emergence_results = get_all_emergence(network, model_calc_params, measure_params, ...
		micro_variable_names, macro_variable_names, pathin, ...
		'measure_params_phiid_ce', measure_params_phiid_ce, ...
		'measure_params_dd', measure_params_dd);

% save emergence_results...
load('/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/12node_kuramoto/01_emergence_results.mat', 'emergence_results');
%}