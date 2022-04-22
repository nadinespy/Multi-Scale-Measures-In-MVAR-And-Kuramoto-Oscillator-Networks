function get_all_kuramoto_log_variables(network, A_vec, beta_vec, all_npoints, pathout_data_sim_time_series)
	
	constant = 100;
	
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
				load([pathout_data_sim_time_series network '_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'pair_sync');
				load([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'mean_pair_sync');
				load([pathout_data_sim_time_series network '_shuffled_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'shuffled_sigma_chi');
				load([pathout_data_sim_time_series network '_shuffled_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'shuffled_phase');
				
				% binarize all variables (add constant to all values 
				% so that there are no negative values)
				log_phase = log(phase+constant);
				log_raw_signal = log(raw_signal+constant);
				log_synchrony = log(synchrony+constant);
				log_pair_sync = log(pair_sync+constant);
				log_mean_pair_sync = log(mean_pair_sync+constant);
				log_sigma_chi = log(sigma_chi+constant);
				log_shuffled_sigma_chi = log(shuffled_sigma_chi+constant);
				log_shuffled_phase = log(shuffled_phase+constant);
				
				% load simulated model with given A and beta
				save([pathout_data_sim_time_series network '_log_phase_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'log_phase');
				save([pathout_data_sim_time_series network '_log_raw_signal_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'log_raw_signal');
				save([pathout_data_sim_time_series network '_log_sigma_chi_' A_str '_' beta_str '_' num2str(npoints)  '.mat'], ...
					'log_sigma_chi');
				save([pathout_data_sim_time_series network '_log_synchrony_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'log_synchrony');
				save([pathout_data_sim_time_series network '_log_mean_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'log_mean_pair_sync');
				save([pathout_data_sim_time_series network '_log_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					'log_pair_sync');
				save([pathout_data_sim_time_series network '_log_shuffled_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'log_shuffled_sigma_chi');
				save([pathout_data_sim_time_series network '_log_shuffled_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'log_shuffled_phase');
				
			end
		end
	end
	
end