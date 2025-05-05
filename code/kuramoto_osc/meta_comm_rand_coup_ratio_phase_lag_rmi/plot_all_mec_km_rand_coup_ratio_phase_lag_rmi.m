%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_all_mec_km_rand_coup_ratio_phase_lag_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In order to run this script, 
% - params_phiid_mmi_rand_coup_ratio_phase_lag_rmi.m
% - mec_km_rand_coup_ratio_phase_lag_rmi_local_dirs.m
% (or any other script with the prefix 'params') must be run first
% (in that order).
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(['/media/nadinespy/NewVolume1/work/phd/projects/mec_mvar_km/', ...
	'mec_simulations/code/kuramoto_osc/' ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi'])

phiid_mmi_measure_in_filename = {'phiid_mmi_ce', 'phiid_mmi_dc', 'phiid_mmi_cd', ...
	'phiid_mmi_uc', 'phiid_mmi_double_red', 'phiid_mmi_syn', 'phiid_mmi_transfer'}; 

phiid_ccs_measure_in_filename = {'phiid_ccs_ce', 'phiid_ccs_dc', ...
	'phiid_ccs_cd', 'phiid_ccs_uc', 'phiid_ccs_double_red', ...
	'phiid_ccs_syn', 'phiid_ccs_transfer'};

for index1 = 1:length(measure_in_simulation)

	if strcmp(measure_in_simulation{index1}, 'phiid_mmi')
		% define measure names as in filenames
		for index2 = 1:length(phiid_mmi_measure_in_filename)
			measure_in_filename = {phiid_mmi_measure_in_filename{index2}};
			plot_mec_km_rand_coup_ratio_phase_lag_rmi();
		end
		
	elseif strcmp(measure_in_simulation{index1}, 'phiid_ccs')
		% define measure names as in filenames
		for index2 = 1:length(phiid_ccs_measure_in_filename)
			measure_in_filename = {phiid_ccs_measure_in_filename{index2}};
			plot_mec_km_rand_coup_ratio_phase_lag_rmi();
		end

	elseif strcmp(measure_in_simulation{index1}, 'dd_ce_co_info') || ...
			strcmp(measure_in_simulation{index1}, 'multi_info') || ...
			strcmp(measure_in_simulation{index1}, 'integrated_info') || ... 
			strcmp(measure_in_simulation{index1}, 'control')		
		% define measure name as in filename
		measure_in_filename = {measure_in_simulation{index1}};
		plot_mec_km_rand_coup_ratio_phase_lag_rmi();
	end
end