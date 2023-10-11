%% directories for MVAR models

mvgc_path	= getenv('MVGC2_PATH');
gpmat_path  = getenv('GPMAT_PATH');
ssdi_path	= getenv('SSDI_PATH');
run(fullfile(mvgc_path,'startup'));

ssdi_path = '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ssdi/';
run(fullfile(ssdi_path,'startup'));

cd '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/common'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/analytical/functions'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/analytical/scripts'
addpath '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/analytical/scripts/measures_random_couplings'
addpath '/media/nadinespy/NewVolume1/work/packages_and_code_repos/ssdi/'

% infodynamics.jar must be (also) located in the same folder as the script that is using it
javaaddpath('/common/infodynamics.jar');

pathout_plots_measures	= ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/results/plots/' num2str(n_nodes) 'node_mvar/analytical/measures_random_couplings/'];
pathout_data_measures	= ['/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/results/analyses/' num2str(n_nodes) 'node_mvar/analytical/measures_random_couplings/'];


