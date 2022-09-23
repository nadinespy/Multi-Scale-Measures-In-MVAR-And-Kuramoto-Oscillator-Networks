%% TO DO



%% create coupling matrices

% {

coupling_matrices = get_coupling_matrices(get_coupling_matrix, model_sim_params);
									
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

% get_kuramoto_variables() is the script to generate all micro and 
% macro variables; modify according to what micro and macro variables you 
% wish to have

% get_all_kuramoto_variables() will generate 'phase', 'raw', 'sync', 
% p_sync', 'mp_sync' and 'chi' as variables

get_variables(network, model_sim_params, coupling_matrices, pathout);

%}

%% get mean of variance of synchronies across time & mean of variance of synchronies across communities

%{

get_km_met_chi(network, npoints, model_params, ...
		pathout_data_sync)
	
%}

%% quantilize micro and macro variables 

% { 

% generated file names consist of 
% [n_oscillators + network + '_quant' + bin + micro/macro variable name + measure_params{1} + measure_params{2} + time_series_length]
get_all_quant_variables(network, measure_params, model_sim_params, ...
	micro_variable_names, macro_variable_names, pathin_sim_time_series);

% discretizing
% get_all_disc_variables(network, model_param1, model_param2, npoints, micro_variable_names, ...
% 	macro_variable_names, bin_number, pathout_data_sim_time_series);

% binarizing
% get_all_kuramoto_bin_variables(network, model_param1, model_param2, npoints, bin_threshold_phase, ...
% 		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
% 		bin_threshold_sigma_chi, pathout_data_sim_time_series);

%}
