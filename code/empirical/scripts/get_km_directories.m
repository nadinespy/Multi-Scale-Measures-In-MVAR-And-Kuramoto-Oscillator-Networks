%% initialize paths for n kuramoto oscillators

cd '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/common'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/empirical/functions'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/empirical/scripts'
addpath '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ReconcilingEmergences-master'
addpath(genpath('/media/nadinespy/NewVolume1/work/packages_and_code_repos/kuramoto'))
addpath(genpath('/media/nadinespy/NewVolume1/work/packages_and_code_repos/mutual_info_kNN'))

% infodynamics.jar must be (also) located in the same folder as the script that is using it
javaaddpath('/common/infodynamics.jar');


% directories for generated data/results
pathout_data = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/' ...
	'results/analyses/' n_oscillators 'node_km/empirical/'];
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/' ...
	'mec_simulations/results/analyses/' n_oscillators 'node_km/empirical/sim_time_series/'];
pathout_data_measures = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/' ...
	'results/analyses/' n_oscillators 'node_km/empirical/measures/'];
pathout_data_metastability = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/' ...
	'mec_simulations/results/analyses/' n_oscillators 'node_km/empirical/metastability/'];
pathout_data_sync = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/' ...
	'results/analyses/' n_oscillators 'node_km/empirical/synchronies/'];

% directories for generated plots
pathout_plots_measures = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/' ...
	'results/plots/' n_oscillators 'node_km/empirical/measures/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/' ...
	'mec_simulations/results/plots/' n_oscillators 'node_km/empirical/distributions/'];
pathout_plots_metastability = ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/' ...
	'mec_simulations/results/plots/' n_oscillators 'node_km/empirical/metastability/'];






