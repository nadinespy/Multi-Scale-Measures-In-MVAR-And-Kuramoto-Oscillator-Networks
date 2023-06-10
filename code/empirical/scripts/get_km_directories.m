%% initialize paths for n kuramoto oscillators

cd '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/empirical/'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/empirical/functions/'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/empirical/scripts/'
addpath '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

% directories for generated data/results
pathout_data = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_oscillators 'node_km/empirical/'];
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_oscillators 'node_km/empirical/sim_time_series/'];
pathout_data_measures = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_oscillators 'node_km/empirical/measures/'];
pathout_data_metastability = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_oscillators 'node_km/empirical/metastability/'];
pathout_data_sync = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/' n_oscillators 'node_km/empirical/synchronies/'];

% directories for generated plots
pathout_plots_measures = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/plots/' n_oscillators 'node_km/empirical/measures/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/plots/' n_oscillators 'node_km/empirical/distributions/'];
pathout_plots_metastability = ['/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/' ...
	'emergence_complexity_simulations/emergence_complexity_simulations/results/plots/' n_oscillators 'node_km/empirical/metastability/'];






