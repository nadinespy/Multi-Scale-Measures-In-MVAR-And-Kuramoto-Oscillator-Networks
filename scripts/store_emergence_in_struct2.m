function [ce_ccs, ce_mmi, ce_practical] = store_emergence_in_struct2(...
			ce_mmi_theta, ce_mmi_cos_theta, ce_mmi_sync, ce_mmi_bin_sync,...
			cd_mmi_theta, cd_mmi_cos_theta, cd_mmi_sync, cd_mmi_bin_sync,...
			dc_mmi_theta, dc_mmi_cos_theta, dc_mmi_sync, dc_mmi_bin_sync, ...
			ce_ccs_theta, ce_ccs_cos_theta, ce_ccs_sync, ce_ccs_bin_sync,...
			cd_ccs_theta, cd_ccs_cos_theta, cd_ccs_sync, cd_ccs_bin_sync,...
			dc_ccs_theta, dc_ccs_cos_theta, dc_ccs_sync, dc_ccs_bin_sync, ...
			ce_pract_sigma_chi_micro_theta, ce_pract_sigma_chi_micro_cos_theta, ce_pract_sigma_chi_micro_sync, ce_pract_sigma_chi_micro_bin_sync,...
			cd_pract_sigma_chi_micro_theta, cd_pract_sigma_chi_micro_cos_theta, cd_pract_sigma_chi_micro_sync, cd_pract_sigma_chi_micro_bin_sync,...
			dc_pract_sigma_chi_micro_theta, dc_pract_sigma_chi_micro_cos_theta, dc_pract_sigma_chi_micro_sync, dc_pract_sigma_chi_micro_bin_sync, ...
			ce_pract_pair_sync_micro_theta, ce_pract_pair_sync_micro_cos_theta, ce_pract_pair_sync_micro_sync, ce_pract_pair_sync_micro_bin_sync,...
			cd_pract_pair_sync_micro_theta, cd_pract_pair_sync_micro_cos_theta, cd_pract_pair_sync_micro_sync, cd_pract_pair_sync_micro_bin_sync,...
			dc_pract_pair_sync_micro_theta, dc_pract_pair_sync_micro_cos_theta, dc_pract_pair_sync_micro_sync, dc_pract_pair_sync_micro_bin_sync)
	
	
	ce_ccs = [];
	ce_mmi = [];
	ce_practical = [];
	
	ce_ccs.ce_ccs_theta = ce_ccs_theta;
	ce_ccs.cd_ccs_theta = cd_ccs_theta;
	ce_ccs.dc_ccs_theta = dc_ccs_theta;
	ce_mmi.ce_mmi_theta = ce_mmi_theta;
	ce_mmi.cd_mmi_theta = cd_mmi_theta;
	ce_mmi.dc_mmi_theta = dc_mmi_theta;
	
	ce_ccs.ce_ccs_cos_theta = ce_ccs_cos_theta;
	ce_ccs.cd_ccs_cos_theta = cd_ccs_cos_theta;
	ce_ccs.dc_ccs_cos_theta = dc_ccs_cos_theta;
	ce_mmi.ce_mmi_cos_theta = ce_mmi_cos_theta;
	ce_mmi.cd_mmi_cos_theta = cd_mmi_cos_theta;
	ce_mmi.dc_mmi_cos_theta = dc_mmi_cos_theta;
	
	ce_ccs.ce_ccs_sync = ce_ccs_sync;
	ce_ccs.cd_ccs_sync = cd_ccs_sync;
	ce_ccs.dc_ccs_sync = dc_ccs_sync;
	ce_mmi.ce_mmi_sync = ce_mmi_sync;
	ce_mmi.cd_mmi_sync = cd_mmi_sync;
	ce_mmi.dc_mmi_sync = dc_mmi_sync;
	
	ce_ccs.ce_ccs_bin_sync = ce_ccs_bin_sync;
	ce_ccs.cd_ccs_bin_sync = cd_ccs_bin_sync;
	ce_ccs.dc_ccs_bin_sync = dc_ccs_bin_sync;
	ce_mmi.ce_mmi_bin_sync = ce_mmi_bin_sync;
	ce_mmi.cd_mmi_bin_sync = cd_mmi_bin_sync;
	ce_mmi.dc_mmi_bin_sync = dc_mmi_bin_sync;
	

	ce_practical.ce_pract_sigma_chi_micro_theta = ce_pract_sigma_chi_micro_theta;
	ce_practical.ce_pract_sigma_chi_micro_cos_theta = ce_pract_sigma_chi_micro_cos_theta;
	ce_practical.ce_pract_sigma_chi_micro_sync = ce_pract_sigma_chi_micro_sync;
	ce_practical.ce_pract_sigma_chi_micro_bin_sync = ce_pract_sigma_chi_micro_bin_sync;

	ce_practical.cd_pract_sigma_chi_micro_theta = cd_pract_sigma_chi_micro_theta;
	ce_practical.cd_pract_sigma_chi_micro_cos_theta = cd_pract_sigma_chi_micro_cos_theta;
	ce_practical.cd_pract_sigma_chi_micro_sync = cd_pract_sigma_chi_micro_sync;
	ce_practical.cd_pract_sigma_chi_micro_bin_sync = cd_pract_sigma_chi_micro_bin_sync;

	ce_practical.dc_pract_sigma_chi_micro_theta = dc_pract_sigma_chi_micro_theta;
	ce_practical.dc_pract_sigma_chi_micro_cos_theta = dc_pract_sigma_chi_micro_cos_theta;
	ce_practical.dc_pract_sigma_chi_micro_sync = dc_pract_sigma_chi_micro_sync;
	ce_practical.dc_pract_sigma_chi_micro_bin_sync = dc_pract_sigma_chi_micro_bin_sync;
	
	
	ce_practical.ce_pract_pair_sync_micro_theta = ce_pract_pair_sync_micro_theta;
	ce_practical.ce_pract_pair_sync_micro_cos_theta = ce_pract_pair_sync_micro_cos_theta;
	ce_practical.ce_pract_pair_sync_micro_sync = ce_pract_pair_sync_micro_sync;
	ce_practical.ce_pract_pair_sync_micro_bin_sync = ce_pract_pair_sync_micro_bin_sync;

	ce_practical.cd_pract_pair_sync_micro_theta = cd_pract_pair_sync_micro_theta;
	ce_practical.cd_pract_pair_sync_micro_cos_theta = cd_pract_pair_sync_micro_cos_theta;
	ce_practical.cd_pract_pair_sync_micro_sync = cd_pract_pair_sync_micro_sync;
	ce_practical.cd_pract_pair_sync_micro_bin_sync = cd_pract_pair_sync_micro_bin_sync;

	ce_practical.dc_pract_pair_sync_micro_theta = dc_pract_pair_sync_micro_theta;
	ce_practical.dc_pract_pair_sync_micro_cos_theta = dc_pract_pair_sync_micro_cos_theta;
	ce_practical.dc_pract_pair_sync_micro_sync = dc_pract_pair_sync_micro_sync;
	ce_practical.dc_pract_pair_sync_micro_bin_sync = dc_pract_pair_sync_micro_bin_sync;

end