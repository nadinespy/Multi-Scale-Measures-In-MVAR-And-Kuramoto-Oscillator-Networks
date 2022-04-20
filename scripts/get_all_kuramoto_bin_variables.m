function get_all_kuramoto_bin_variables(network, A_vec, beta_vec, all_npoints, bin_threshold_phase, ...
		bin_threshold_raw_signal, bin_threshold_sync, bin_threshold_pair_sync, ...
		bin_threshold_sigma_chi, pathout_data_sim_time_series);
	
	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		rng(1);
		for i = 1:length(A_vec);
			A_str = param2str(A_vec(i));
			
			for j = 1:length(beta_vec)
				beta_str = param2str(beta_vec(j));
				
				% load simulated model with given A and beta:
				% micro variables - phases, raw signal, synchronies;
				% macro variables for practical measures for causal emergence -
				% variance of synchronies & global average pairwise synchrony
				% between communities
				load([pathout_data_sim_time_series network '_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'phase');
				load([pathout_data_sim_time_series network '_synchrony_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'synchrony');
				load([pathout_data_sim_time_series network '_raw_signal_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'raw_signal');
				load([pathout_data_sim_time_series network '_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'sigma_chi');
				load([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'grand_mean_pair_sync');
				
				% binarize all variables
				bin_phase = get_binarized_variables(phase, bin_threshold_phase);
				bin_raw_signal = get_binarized_variables(raw_signal, bin_threshold_raw_signal);
				bin_synchrony = get_binarized_variables(synchrony, bin_threshold_sync);
				bin_grand_mean_pair_sync = get_binarized_variables(grand_mean_pair_sync, bin_threshold_pair_sync);
				bin_sigma_chi = get_binarized_variables(sigma_chi, bin_threshold_sigma_chi);
				
				% load simulated model with given A and beta
				save([pathout_data_sim_time_series network '_bin_phase_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'bin_phase');
				save([pathout_data_sim_time_series network '_bin_raw_signal_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'bin_raw_signal');
				save([pathout_data_sim_time_series network '_bin_sigma_chi_' A_str '_' beta_str '_' num2str(npoints)  '.mat'], ...
					'bin_sigma_chi');
				save([pathout_data_sim_time_series network '_bin_synchrony_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'bin_synchrony');
				save([pathout_data_sim_time_series network '_bin_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'bin_grand_mean_pair_sync');
				
			end
		end
	end
	
end