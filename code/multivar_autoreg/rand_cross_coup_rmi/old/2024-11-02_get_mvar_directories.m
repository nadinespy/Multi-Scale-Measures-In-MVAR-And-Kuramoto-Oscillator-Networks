

%% directories for MVAR models

mvgc_path	= getenv('MVGC2_PATH');
gpmat_path	= getenv('GPMAT_PATH');
ssdi_path	= getenv('SSDI_PATH');
ssdi_path	= '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ssdi/';

run(fullfile(mvgc_path,'startup'));
run(fullfile(ssdi_path,'startup'));

basedir = ['/media/nadinespy/NewVolume1/work/phd/projects/' ...
	'mec_experiments/mec_simulations/'];

addpath([basedir 'code']);
addpath(basedir 'code/common']);
addpath([basedir 'code/analytical/functions']);
addpath([basedir 'code/analytical/scripts']);
addpath([basedir 'code/analytical/scripts/mec_random_cross_coups_rmi'];)
addpath '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ssdi/'

cd '/media/nadinespy/NewVolume1/work/phd/projects/mec_experiments/mec_simulations/code'

% infodynamics.jar must be (also) located in the same folder as the script that is using it
javaaddpath('/common/infodynamics.jar');

%pathout_plots_measures	= ['/media/nadinespy/NewVolume1/work/phd/projects/mec_experiments/mec_simulations/results/plots/' ...
%	num2str(n_nodes) 'node_mvar/analytical/measures_random_couplings/'];
pathout_plots_measures	= ['/media/nadinespy/UUI/Nadine/current_results_plots_phd/multivar_autoreg/' ...
	'8mvar_random_cross_coup/analytical/mec_random_cross_coups_rmi'];

%pathout_data_measures	= ['/media/nadinespy/NewVolume1/work/phd/projects/mec_experiments/mec_simulations/results/analyses/' ...
%	num2str(n_nodes) 'node_mvar/analytical/measures_random_couplings/'];
pathout_data_measures	= ['/media/nadinespy/UUI/Nadine/current_results_phd/multivar_autoreg/' ...
	'8mvar_random_cross_coup/analytical/mec_random_cross_coups_rmi'];

