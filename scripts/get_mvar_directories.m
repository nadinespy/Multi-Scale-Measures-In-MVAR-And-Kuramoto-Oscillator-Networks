%% initialize paths for n nodes

cd '/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/EmergenceComplexityMeasuresComparison_Matlab/scripts'
addpath '/media/nadinespy/NewVolume/work/toolboxes_matlab'
addpath '/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/EmergenceComplexityMeasuresComparison_Matlab/scripts/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

% directories for generated data/results
pathout_data = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/' n_nodes 'node_mvar/'];
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume1//work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/' n_nodes 'node_mvar/sim_time_series/'];
pathout_data_measures = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/' n_nodes 'node_mvar/measures/'];
pathout_data_mean_corr = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/' n_nodes 'node_mvar/average_corr/'];
pathout_data_mean_cov = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/' n_nodes 'node_mvar/average_cov/'];

% directories for generated plots
pathout_plots_measures = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/' n_nodes 'node_mvar/measures/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/' n_nodes 'node_mvar/distributions/'];
pathout_plots_mean_corr = ['/media/nadinespy/NewVolume1/work/PhD/my_projects/EmergenceComplexityMeasuresComparisonSimulations/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/' n_nodes 'node_mvar/mean_corr/'];

