%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_all_mec_mvar_rand_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In order to run this script, first run in the given order
% - params_[some_measure_from_mvar_sims]_rand_cross_coup_rmi.m
% - mec_mvar_rand_cross_coup_rmi_local_dirs.m
% - mec_mvar_rand_cross_coup_rmi.m (if results are not generated yet)
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(['/media/nadinespy/NewVolume1/work/phd/projects/mec_mvar_km/', ...
	'mec_simulations/code/multivar_autoreg/', 'rand_cross_coup_rmi']);

phiid_mmi_measure_in_filename = {'phiid_mmi_ce', 'phiid_mmi_dc', 'phiid_mmi_cd', ...
	'phiid_mmi_uc', 'phiid_mmi_double_red', 'phiid_mmi_syn', 'phiid_mmi_transfer'}; 

phiid_ccs_measure_in_filename = {'phiid_ccs_ce', 'phiid_ccs_dc', ...
	'phiid_ccs_cd', 'phiid_ccs_uc', 'phiid_ccs_double_red', ...
	'phiid_ccs_syn', 'phiid_ccs_transfer'};

for index1 = 1:length(measure_in_simulation)

	if strcmp(measure_in_simulation{index1}, 'phiid_mmi')
		% define measure name as in filename
		for index2 = 1:length(phiid_mmi_measure_in_filename)
			measure_in_filename = {phiid_mmi_measure_in_filename{index2}};
			plot_mec_mvar_rand_cross_coup_rmi();
		end
		
	elseif strcmp(measure_in_simulation{index1}, 'phiid_ccs')
		% define measure name as in filename
		for index2 = 1:length(phiid_ccs_measure_in_filename)
			measure_in_filename = {phiid_ccs_measure_in_filename{index2}};
			plot_mec_mvar_rand_cross_coup_rmi();
		end

	elseif strcmp(measure_in_simulation{index1}, 'dd_ce_co_info') || ...
			strcmp(measure_in_simulation{index1}, 'integrated_info') || ... 
			strcmp(measure_in_simulation{index1}, 'control')		
		% define measure name as in filename
		measure_in_filename = {measure_in_simulation{index1}};
		plot_mec_mvar_rand_cross_coup_rmi();
	end

end
