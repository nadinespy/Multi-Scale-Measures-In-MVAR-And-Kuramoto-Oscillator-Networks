%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dirs_and_params_for_local_machine.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Before running this script, run 'params_[measure]_rand_cross_coup_rmi.m' 
% ([measure] as used in parameter specification scripts which start with
% the prefix 'params'.
% 
% This script must be run in order to run 
% - mec_mvar_rand_cross_coup_rmi.m
% - mec_mvar_rand_cross_coup_rmi_plotting.m
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measure (as in the filenames - not in the Matlab structs)
measure_in_filename = 'phiid_ccs_dc';

% directories for MVAR models
basedir = ['/media/nadinespy/NewVolume1/work/phd/projects/' ...
	'mec_experiments/mec_simulations'];

mvgc_path	= getenv('MVGC2_PATH');
gpmat_path	= getenv('GPMAT_PATH');
ssdi_path	= getenv('SSDI_PATH');
ssdi_path	= ['/media/nadinespy/NewVolume1/work/' ...
	'packages_and_code_repos/ssdi/'];

run(fullfile(mvgc_path,'startup'));
run(fullfile(ssdi_path,'startup'));

addpath(fullfile('/media/nadinespy/NewVolume1/work', 'packages_and_code_repos', 'ReconcilingEmergences-master'));
addpath([basedir '/code']);
addpath([basedir '/code/common']);
addpath([basedir '/code/multivar_autoreg/rand_cross_coup_rmi/functions']);
addpath([basedir '/code/multivar_autoreg/rand_cross_coup_rmi/scripts']);

cd '/media/nadinespy/NewVolume1/work/phd/projects/mec_experiments/mec_simulations/code/'

% infodynamics.jar must be (also) located in the same folder as the script that is using it
javaaddpath('/common/infodynamics.jar');

pathout_plots_measures	= [basedir '/results/plots/multivar_autoreg/rand_cross_coup_rmi/' ...
	model_specific_path '/'];
pathout_data_measures	= [basedir '/results/analyses/multivar_autoreg/rand_cross_coup_rmi/' ...
	model_specific_path '/'];

% for plotting
all_cross_coup_mag = 2.^linspace(cross_coup_range(1),cross_coup_range(2), n_cross_coup); % norm values (exponential scale)
all_rmi   = 2.^linspace(rmi_range(1), rmi_range(2), n_rmi); % rmi values (exponential scale)




