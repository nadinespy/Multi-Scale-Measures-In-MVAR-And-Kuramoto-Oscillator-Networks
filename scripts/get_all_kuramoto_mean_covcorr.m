function get_all_kuramoto_mean_covcorr(network, A_vec, beta_vec, all_npoints, ...
		pathout_data_sim_time_series, ...
		pathout_data_mean_corr, ...
		pathout_data_mean_cov);
	
	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		for i = 1:length(A_vec);
			A_str = param2str(A_vec(i));
			
			for j = 1:length(beta_vec)
				beta = beta_vec(j);
				beta_str = param2str(beta);
				
				% load simulated model with given A and beta:
				% micro variables - phases, raw signal, synchronies
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
				
				% MEAN COVARIANCES/CORRELATIONS (calculate mean MI as opposed to correlation?)
				
				% between micro variables
				[mean_cov_phase(i,j), mean_corr_phase(i,j)] = get_mean_covcorr_one_matrix(phase');
				[mean_cov_raw_signal(i,j), mean_corr_raw_signal(i,j)] = get_mean_covcorr_one_matrix(raw_signal');
				[mean_cov_synchrony(i,j), mean_corr_synchrony(i,j)] = get_mean_covcorr_one_matrix(synchrony');
				
				% between micro and macro variables
				[mean_cov_phase_sigma_chi(i,j), mean_corr_phase_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(phase', sigma_chi');
				[mean_cov_phase_pair_sync(i,j), mean_corr_phase_pair_sync(i,j)] = get_mean_covcorr_two_matrices(phase', grand_mean_pair_sync');
				
				[mean_cov_raw_signal_sigma_chi(i,j), mean_corr_raw_signal_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(raw_signal', sigma_chi');
				[mean_cov_raw_signal_pair_sync(i,j), mean_corr_raw_signal_pair_sync(i,j)] = get_mean_covcorr_two_matrices(raw_signal', grand_mean_pair_sync');
				
				[mean_cov_synchrony_sigma_chi(i,j), mean_corr_synchrony_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(synchrony', sigma_chi');
				[mean_cov_synchrony_pair_sync(i,j), mean_corr_synchrony_pair_sync(i,j)] = get_mean_covcorr_two_matrices(synchrony', grand_mean_pair_sync');
				
				
			end
		end
		
		% save covariances and correlations for all values of A and all values of beta; saved filenames consist of
		% network name + '_mean_' + 'corr_' or 'cov_' + micro variable name (+ macro variable name) + number of datapoints
		
		% micro variables
		save([pathout_data_mean_corr network '_mean_corr_phase_' num2str(npoints) '.mat'], ...
			'mean_corr_phase');
		save([pathout_data_mean_cov network '_mean_cov_phase_' num2str(npoints) '.mat'], ...
			'mean_cov_phase');
		
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal');
		
		save([pathout_data_mean_corr network '_mean_corr_synchrony_' num2str(npoints) '.mat'], ...
			'mean_corr_synchrony');
		save([pathout_data_mean_cov network '_mean_cov_synchrony_' num2str(npoints) '.mat'], ...
			'mean_cov_synchrony');
		
		% micro and macro variables
		
		% phase
		save([pathout_data_mean_corr network '_mean_corr_phase_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_phase_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_phase_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_phase_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_pair_sync');
		
		% raw signal
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_pair_sync');
		
		% synchrony
		save([pathout_data_mean_corr network '_mean_corr_synchrony_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_synchrony_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_synchrony_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_synchrony_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_synchrony_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_synchrony_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_synchrony_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_synchrony_pair_sync');
		
	end

end 