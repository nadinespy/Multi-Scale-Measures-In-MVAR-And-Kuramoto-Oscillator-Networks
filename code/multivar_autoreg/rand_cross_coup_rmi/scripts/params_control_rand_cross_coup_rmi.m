% parameters for 'mec_mvar_rand_cross_coup_rmi'

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
rmi_range                       = [-1.5, 1.5];
n_rmi                           = 100;
seed                            = 1;

% --------------------------------------------
% measure-specific parameters
% --------------------------------------------

% list of integrated info measures (in this case,
% control measures)
list_integrated_info_measures = {'TimeDelayedMutualInfo', 'AverageCorrelation'};

% first part of filename
results_filename        = 'control';

% --------------------------------------------
% model-specific parameters
% --------------------------------------------

% second part of filename (must include # of nodes and lag)
model_specific_filename = '8mvar_lag1_rand_cross_coup_rmi';

% model-specific directory (must include # of nodes and lag)
model_specific_path     = '8mvar_lag1';
