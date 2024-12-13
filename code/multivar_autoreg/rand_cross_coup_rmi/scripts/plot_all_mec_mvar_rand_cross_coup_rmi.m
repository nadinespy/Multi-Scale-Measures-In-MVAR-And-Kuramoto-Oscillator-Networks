%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_all_mec_mvar_rand_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In order to run this script, 
% - dirs_and_params_for_local_machine.m
% - params_phiid_mmi_rand_cross_coup_rmi.m 
% (or any other script with the prefix 'params') must be run.
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define array of measure names as used in simulations
measure_in_simulation = {'phiid_mmi', 'phiid_ccs', 'dd_ce_co_info', 'multi_info', ...
	'integrated_info', 'control'};

for index1 = 1:length(measure_in_simulation)

	if strcmp(measure_in_simulation{index1}, 'phiid_mmi')
		% define measure names as in filenames
		measure_in_filename = {'phiid_mmi_ce', 'phiid_mmi_dc', ...
			'phiid_mmi_cd', 'phiid_mmi_uc', 'phiid_mmi_double_red', ...
			'phiid_mmi_syn', 'phiid_mmi_transfer'};
		
	elseif strcmp(measure_in_simulation{index1}, 'phiid_ccs')
		% define measure names as in filenames
		measure_in_filename = {'phiid_ccs_ce', 'phiid_ccs_dc', ...
			'phiid_ccs_cd', 'phiid_ccs_uc', 'phiid_ccs_double_red', ...
			'phiid_ccs_syn', 'phiid_ccs_transfer'};

	elseif strcmp(measure_in_simulation{index1}, 'dd_ce_co_info') || ...
			strcmp(measure_in_simulation{index1}, 'multi_info') || ...
			strcmp(measure_in_simulation{index1}, 'integrated_info') || ... 
			strcmp(measure_in_simulation{index1}, 'control')		
		% define measure name as in filename
		measure_in_filename = {measure_in_simulation{index1}};

	for index2 = 1:length(measure_in_filename)
		plot_mec_mvar_rand_cross_coup_rmi();
	end
end