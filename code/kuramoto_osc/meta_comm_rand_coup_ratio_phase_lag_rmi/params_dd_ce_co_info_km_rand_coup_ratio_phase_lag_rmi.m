%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% params_dd_ce_co_info_km_rand_coup_ratio_phase_lag_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This scripts defines parameters for 'mec_km_rand_coup_ratio_phase_lag_rmi'.
% 
% This is the first script to run in order to then run
% - mec_km_rand_coup_ratio_phase_lag_rmi_local_dirs.m
% - mec_km_rand_coup_ratio_phase_lag_rmi.m
% - plot_all_mec_km_rand_coup_ratio_phase_lag_rmi.m
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model order for DD calculation
ar_model_order			= 1;

% network topology
n_comms				= 4;			% number of communities
n_inter_comm_coups		= 8;			% number of oscillators per community
n_nodes_per_comm			= 8;			% number of inter-group couplings per node
n_nodes				= n_comms * n_nodes_per_comm; % total number of oscillators

% oscillator parameters
osc_freq_mean			= 1;			% oscillator frequencies mean
osc_freq_rel_dev			= 0;			% oscillator frequencies relative deviation

phase_lag_range			= [0, 0.8];		
n_phase_lag				= 100;

% coupling strength parameters
n_samples_coup			= 1;			% number of samples of coupling matrices: if > 1, 
								% coup_ratio_mean_rel_dev needs to be > 0!
coup_ratio_mean_range		= [0, 1];
n_coup_ratio_mean			= 100;		% number of coupling ratios
coup_ratio_mean_rel_dev		= 0;			% connectivity parameter relative deviation: if > 0,
								% n_samples_coup_ratio_mean needs to be > 1!
global_coup				= 1;			% global connectivity factor

% noise correlation parameters
n_samples_noise_corr		= 1;			% number of samples of noise correlation matrices 
rmi_range				= [0, 0];
n_rmi					= 1;
seed					= 1;

% time parameters
n_sim_points			= 100;		% equilibriated simulation time
n_equil_points			= 100;		% equilibriation time
dt_sim_points			= 0.01;		% integration time increment
sim_mode				= 'Euler';		% simulation mode: 'Euler' or 'RK4'

% noise correlation parameters
noise_mag				= 0.1;		% oscillator input noise magnitude
adjust_rmi				= true;		% oscillator input noise multi-information 
								% adjusted for system size?

% --------------------------------------------
% measure-specific parameters
% --------------------------------------------

% which macro variables & number of dimensions 

dim_reduction         	= {'comm_sync'}; % 'full_system_sync' ,'chim_index'};
m_dim_ce_dd_co_info	= n_comms;

% VOU-DD sampling parameters
n_local_gc			= 1000;		% local Granger-Causality sample size
emp_sample			= false;		% empirical sample (else uniform over 
							% oscillator network state space)?
stabilise_care		= sqrt(eps);	% stabilising perturbation for CARE


% first part of filename
results_filename        = 'dd_ce_co_info';

% --------------------------------------------
% model-specific parameters
% --------------------------------------------

% second part of filename (must include # of nodes and lags)
model_specific_filename = [num2str(n_nodes) 'km_lag' num2str(ar_model_order) ...
	'_rand_coup_ratio_phase_lag_rmi'];

% model-specific directory (must include # of nodes and lags)
model_specific_path     = [num2str(n_nodes) 'km_lag' num2str(ar_model_order)];
