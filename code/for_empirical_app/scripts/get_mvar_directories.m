%% initialize paths for n nodes

cd '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/empirical/'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/'
addpath '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

% directories for generated data/results
pathout_data = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_nodes 'node_mvar/empirical/'];
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_nodes 'node_mvar/empirical/sim_time_series/'];
pathout_data_measures = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_nodes 'node_mvar/empirical/measures/'];

% directories for generated plots
pathout_plots_measures = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/plots/' n_nodes 'node_mvar/empirical/measures/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/plots/' n_nodes 'node_mvar/empirical/distributions/'];

