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
				load([pathout_data_sim_time_series network '_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'pair_sync');
				load([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'mean_pair_sync');
				load([pathout_data_sim_time_series network '_shuffled_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'shuffled_sigma_chi');
				load([pathout_data_sim_time_series network '_shuffled_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'shuffled_phase');
				
				% MEAN COVARIANCES/CORRELATIONS (calculate mean MI as opposed to correlation?)
				
				% between micro variables
				[mean_cov_phase(i,j), mean_corr_phase(i,j)] = get_mean_covcorr_one_matrix(phase');
				[mean_cov_raw_signal(i,j), mean_corr_raw_signal(i,j)] = get_mean_covcorr_one_matrix(raw_signal');
				[mean_cov_synchrony(i,j), mean_corr_synchrony(i,j)] = get_mean_covcorr_one_matrix(synchrony');
				[mean_cov_pair_sync(i,j), mean_corr_pair_sync(i,j)] = get_mean_covcorr_one_matrix(pair_sync');
				[mean_cov_shuffled_phase(i,j), mean_corr_shuffled_phase(i,j)] = get_mean_covcorr_one_matrix(shuffled_phase');
				
				% between micro and macro variables
				[mean_cov_phase_sigma_chi(i,j), mean_corr_phase_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(phase', sigma_chi');
				[mean_cov_phase_mean_pair_sync(i,j), mean_corr_phase_mean_pair_sync(i,j)] = get_mean_covcorr_two_matrices(phase', mean_pair_sync');
				[mean_cov_phase_synchrony(i,j), mean_corr_phase_synchrony(i,j)] = get_mean_covcorr_two_matrices(phase', synchrony');
				[mean_cov_phase_pair_sync(i,j), mean_corr_phase_pair_sync(i,j)] = get_mean_covcorr_two_matrices(phase', pair_sync');
				
				[mean_cov_raw_signal_sigma_chi(i,j), mean_corr_raw_signal_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(raw_signal', sigma_chi');
				[mean_cov_raw_signal_mean_pair_sync(i,j), mean_corr_raw_signal_mean_pair_sync(i,j)] = get_mean_covcorr_two_matrices(raw_signal', mean_pair_sync');
				[mean_cov_raw_signal_synchrony(i,j), mean_corr_raw_signal_synchrony(i,j)] = get_mean_covcorr_two_matrices(raw_signal', synchrony');
				[mean_cov_raw_signal_pair_sync(i,j), mean_corr_raw_signal_pair_sync(i,j)] = get_mean_covcorr_two_matrices(raw_signal', pair_sync');
				
				[mean_cov_synchrony_sigma_chi(i,j), mean_corr_synchrony_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(synchrony', sigma_chi');
				[mean_cov_synchrony_mean_pair_sync(i,j), mean_corr_synchrony_mean_pair_sync(i,j)] = get_mean_covcorr_two_matrices(synchrony', mean_pair_sync');
				
				[mean_cov_pair_sync_sigma_chi(i,j), mean_corr_pair_sync_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(pair_sync', sigma_chi');
				[mean_cov_pair_sync_mean_pair_sync(i,j), mean_corr_pair_sync_mean_pair_sync(i,j)] = get_mean_covcorr_two_matrices(pair_sync', mean_pair_sync');
				
				% between shuffled phase and macro variables
				[mean_cov_shuffled_phase_sigma_chi(i,j), mean_corr_shuffled_phase_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(shuffled_phase', sigma_chi');
				[mean_cov_shuffled_phase_mean_pair_sync(i,j), mean_corr_shuffled_phase_mean_pair_sync(i,j)] = get_mean_covcorr_two_matrices(shuffled_phase', mean_pair_sync');
				[mean_cov_shuffled_phase_synchrony(i,j), mean_corr_shuffled_phase_synchrony(i,j)] = get_mean_covcorr_two_matrices(shuffled_phase', synchrony');
				[mean_cov_shuffled_phase_pair_sync(i,j), mean_corr_shuffled_phase_pair_sync(i,j)] = get_mean_covcorr_two_matrices(shuffled_phase', pair_sync');
				
				% between micro variables and shuffled macro variable
				[mean_cov_phase_shuffled_sigma_chi(i,j), mean_corr_phase_shuffled_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(phase', shuffled_sigma_chi');
				[mean_cov_raw_signal_shuffled_sigma_chi(i,j), mean_corr_raw_signal_shuffled_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(raw_signal', shuffled_sigma_chi');
				[mean_cov_synchrony_shuffled_sigma_chi(i,j), mean_corr_synchrony_shuffled_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(synchrony', shuffled_sigma_chi');
				[mean_cov_pair_sync_shuffled_sigma_chi(i,j), mean_corr_pair_sync_shuffled_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(pair_sync', shuffled_sigma_chi');
				
				% between shuffled micro and shuffled macro variable
				[mean_cov_shuffled_phase_shuffled_sigma_chi(i,j), mean_corr_shuffled_phase_shuffled_sigma_chi(i,j)] = get_mean_covcorr_two_matrices(shuffled_phase', shuffled_sigma_chi');
				
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
		
		save([pathout_data_mean_corr network 'mean_corr_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_pair_sync');
		save([pathout_data_mean_cov network 'mean_cov_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_pair_sync');
		
		% micro and macro variables
		
		% phase + macro variables
		save([pathout_data_mean_corr network '_mean_corr_phase_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_phase_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_phase_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_mean_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_phase_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_mean_pair_sync');
		
		save([pathout_data_mean_corr network '_mean_corr_phase_synchrony_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_synchrony');
		save([pathout_data_mean_cov network '_mean_cov_phase_synchrony_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_synchrony');

		save([pathout_data_mean_corr network '_mean_corr_phase_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_phase_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_pair_sync');
		
		% raw signal + macro variables
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_mean_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_mean_pair_sync');
		
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_synchrony_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_synchrony');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_synchrony_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_synchrony');

		save([pathout_data_mean_corr network '_mean_corr_raw_signal_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_pair_sync');
		
		% synchrony + macro variables
		save([pathout_data_mean_corr network '_mean_corr_synchrony_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_synchrony_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_synchrony_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_synchrony_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_synchrony_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_synchrony_mean_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_synchrony_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_synchrony_mean_pair_sync');
		
		% pairwise synchrony + macro variables
		save([pathout_data_mean_corr network '_mean_corr_pair_sync_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_pair_sync_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_pair_sync_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_pair_sync_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_pair_sync_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_pair_sync_mean_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_pair_sync_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_pair_sync_mean_pair_sync');
		
		% shuffled phase + macro variables
		save([pathout_data_mean_corr network '_mean_corr_shuffled_phase_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_shuffled_phase_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_shuffled_phase_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_shuffled_phase_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_shuffled_phase_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_shuffled_phase_mean_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_phase_mean_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_shuffled_phase_mean_pair_sync');
		
		save([pathout_data_mean_corr network '_mean_corr_shuffled_phase_synchrony_' num2str(npoints) '.mat'], ...
			'mean_corr_shuffled_phase_synchrony');
		save([pathout_data_mean_cov network '_mean_cov_phase_synchrony_' num2str(npoints) '.mat'], ...
			'mean_cov_shuffled_phase_synchrony');
		
		save([pathout_data_mean_corr network '_mean_corr_shuffled_phase_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_corr_shuffled_phase_pair_sync');
		save([pathout_data_mean_cov network '_mean_cov_phase_pair_sync_' num2str(npoints) '.mat'], ...
			'mean_cov_shuffled_phase_pair_sync');
		
		% micro variables + shuffled sigma chi
		save([pathout_data_mean_corr network '_mean_corr_phase_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_phase_shuffled_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_phase_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_phase_shuffled_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_raw_signal_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_raw_signal_shuffled_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_raw_signal_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_raw_signal_shuffled_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_synchrony_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_synchrony_shuffled_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_synchrony_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_synchrony_shuffled_sigma_chi');
		
		save([pathout_data_mean_corr network '_mean_corr_pair_sync_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_pair_sync_shuffled_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_pair_sync_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_pair_sync_shuffled_sigma_chi');
		
		% shuffled phase + shuffled sigma chi
		save([pathout_data_mean_corr network '_mean_corr_shuffled_phase_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_corr_shuffled_phase_shuffled_sigma_chi');
		save([pathout_data_mean_cov network '_mean_cov_shuffled_phase_shuffled_sigma_chi_' num2str(npoints) '.mat'], ...
			'mean_cov_shuffled_phase_shuffled_sigma_chi');
		
	end

end 