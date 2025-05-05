%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% params_phiid_ccs_mvar_rand_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This scripts defines parameters for 'mec_mvar_rand_cross_coup_rmi'.
% 
% This is the first script to run in order to then run
% - mec_mvar_rand_cross_coup_rmi_local_dirs.m
% - mec_mvar_rand_cross_coup_rmi.m
% - plot_all_mec_mvar_rand_cross_coup_rmi.m
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------
% same for all measures
% --------------------------------------------
n_nodes                         = 8;
ar_model_order                  = 1;
spectral_radius                 = 0.9;
n_samples_cross_coup            = 25;
n_samples_noise_corr            = 25;
cross_coup_range                = [-5, 5];
n_cross_coup                    = 100;
rmi_range                       = [-2, 2];
n_rmi                           = 100;
seed                            = 1;

% --------------------------------------------
% measure-specific parameters
% --------------------------------------------

% first part of filename
results_filename		= 'phiid_ccs';
n_datapoints_phiid_ccs		= 4000;

% --------------------------------------------
% model-specific parameters
% --------------------------------------------

% second part of filename 
% (must include # of nodes and lag)
model_specific_filename = [num2str(n_nodes) ...
	'mvar_lag' num2str(ar_model_order) ...
	'_rand_cross_coup_rmi'];

% model-specific directory 
% (must include # of nodes and lags)
model_specific_path     = [num2str(n_nodes) ...
	'mvar_lag' num2str(ar_model_order)];

