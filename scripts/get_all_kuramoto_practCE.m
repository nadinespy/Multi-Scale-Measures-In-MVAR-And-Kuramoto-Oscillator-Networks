function get_all_kuramoto_practCE(network, A_vec, beta_vec, coupling_matrices, taus, all_npoints, method, ...
		pathout_data_sim_time_series, pathout_data_pract_ce)

	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		for z = 1:length(taus);
			tau = taus(z);
			
			rng(1);
			for i = 1:size(coupling_matrices, 3);
				A_str = param2str(A_vec(i));
				
				for j = 1:length(beta_vec)
					beta_str = param2str(beta_vec(j));
					
					fprintf('get_all_kuramoto_practCE - loop indices: npoints: %d, tau: %d, A: %d, beta: %d\n', q, z, i, j);
					
					% load simulated model with given A and beta
					load([pathout_data_sim_time_series network '_bin_phase_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'bin_phase');
					load([pathout_data_sim_time_series network '_bin_raw_signal_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'bin_raw_signal');
					load([pathout_data_sim_time_series network '_bin_sigma_chi_' A_str '_' beta_str '_' num2str(npoints)  '.mat'], ...
						'bin_sigma_chi');
					load([pathout_data_sim_time_series network '_bin_synchrony_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'bin_synchrony');
					load([pathout_data_sim_time_series network '_bin_mean_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'bin_mean_pair_sync');
					load([pathout_data_sim_time_series network '_bin_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'bin_pair_sync');
					load([pathout_data_sim_time_series network '_bin_shuffled_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
						'bin_shuffled_sigma_chi');
					load([pathout_data_sim_time_series network '_bin_shuffled_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
						'bin_shuffled_phase');

					
					% micro variables: binarized phase, raw signal, synchrony, pairwise synchrony
					% macro variables: binarized chimera-index (sigma_chi), global average pairwise synchrony
					
					% ---------------------------------------------------------------------------------------------------------------------------------------
					% practical causal emergence
					% ---------------------------------------------------------------------------------------------------------------------------------------
					
					% store micro and macro variables in two structs (variables are transposed, so have time-points in rows)
					macro_variables.bin_sigma_chi = bin_sigma_chi';
					macro_variables.bin_mean_pair_sync = bin_mean_pair_sync';
					macro_variables.bin_shuffled_sigma_chi = bin_shuffled_sigma_chi';

					micro_variables.bin_phase = bin_phase';
					micro_variables.bin_raw_signal = bin_raw_signal';
					micro_variables.bin_synchrony = bin_synchrony';
					micro_variables.bin_pair_sync = bin_pair_sync';
					micro_variables.bin_shuffled_phase = bin_shuffled_phase';
					
					% calculate practical CE
					
					% get practical CE, DC, & CE for all combinations of micro and macro variables
					practCE = get_practCE(micro_variables, macro_variables, tau, method);
					
					% extract practical CE for different combinations of micro and macro variables, and store it in arrays for 
					% different parameter combinations of beta and A
					pract_ce_bin_phase_bin_sigma_chi(i,j) = practCE.bin_phase_bin_sigma_chi.pract_ce;
					pract_ce_bin_phase_bin_mean_pair_sync(i,j) = practCE.bin_phase_bin_mean_pair_sync.pract_ce;
					pract_ce_bin_raw_signal_bin_sigma_chi(i,j) = practCE.bin_raw_signal_bin_sigma_chi.pract_ce;
					pract_ce_bin_raw_signal_bin_mean_pair_sync(i,j) = practCE.bin_raw_signal_bin_mean_pair_sync.pract_ce;
					pract_ce_bin_synchrony_bin_sigma_chi(i,j) = practCE.bin_synchrony_bin_sigma_chi.pract_ce;
					pract_ce_bin_synchrony_bin_mean_pair_sync(i,j) = practCE.bin_synchrony_bin_mean_pair_sync.pract_ce;
					pract_ce_bin_pair_sync_bin_sigma_chi(i,j) = practCE.bin_pair_sync_bin_sigma_chi.pract_ce;
					pract_ce_bin_pair_sync_bin_mean_pair_sync(i,j) = practCE.bin_pair_sync_bin_mean_pair_sync.pract_ce;
					
					% shuffled micro and macro variables
					pract_ce_bin_phase_bin_shuffled_sigma_chi(i,j) = practCE.bin_phase_bin_shuffled_sigma_chi.pract_ce;
					pract_ce_bin_raw_signal_bin_shuffled_sigma_chi(i,j) = practCE.bin_raw_signal_bin_shuffled_sigma_chi.pract_ce;
					pract_ce_bin_synchrony_bin_shuffled_sigma_chi(i,j) = practCE.bin_synchrony_bin_shuffled_sigma_chi.pract_ce;
					pract_ce_bin_pair_sync_bin_shuffled_sigma_chi(i,j) = practCE.bin_pair_sync_bin_shuffled_sigma_chi.pract_ce;
					
					pract_ce_bin_shuffled_phase_bin_sigma_chi(i,j) = practCE.bin_shuffled_phase_bin_sigma_chi.pract_ce;
					pract_ce_bin_shuffled_phase_bin_mean_pair_sync(i,j) = practCE.bin_shuffled_phase_bin_mean_pair_sync.pract_ce;
					
					pract_ce_bin_shuffled_phase_bin_shuffled_sigma_chi(i,j) = practCE.bin_phase_bin_shuffled_sigma_chi.pract_ce;

					% extract practical DC, and store it in arrays for different combinations of micro and macro variables
					pract_dc_bin_phase_bin_sigma_chi(i,j) = practCE.bin_phase_bin_sigma_chi.pract_dc;
					pract_dc_bin_phase_bin_mean_pair_sync(i,j) = practCE.bin_phase_bin_mean_pair_sync.pract_dc;
					pract_dc_bin_raw_signal_bin_sigma_chi(i,j) = practCE.bin_raw_signal_bin_sigma_chi.pract_dc;
					pract_dc_bin_raw_signal_bin_mean_pair_sync(i,j) = practCE.bin_raw_signal_bin_mean_pair_sync.pract_dc;
					pract_dc_bin_synchrony_bin_sigma_chi(i,j) = practCE.bin_synchrony_bin_sigma_chi.pract_dc;
					pract_dc_bin_synchrony_bin_mean_pair_sync(i,j) = practCE.bin_synchrony_bin_mean_pair_sync.pract_dc;
					pract_dc_bin_pair_sync_bin_sigma_chi(i,j) = practCE.bin_pair_sync_bin_sigma_chi.pract_dc;
					pract_dc_bin_pair_sync_bin_mean_pair_sync(i,j) = practCE.bin_pair_sync_bin_mean_pair_sync.pract_dc;
					
					% shuffled micro and macro variables
					pract_dc_bin_phase_bin_shuffled_sigma_chi(i,j) = practCE.bin_phase_bin_shuffled_sigma_chi.pract_dc;
					pract_dc_bin_raw_signal_bin_shuffled_sigma_chi(i,j) = practCE.bin_raw_signal_bin_shuffled_sigma_chi.pract_dc;
					pract_dc_bin_synchrony_bin_shuffled_sigma_chi(i,j) = practCE.bin_synchrony_bin_shuffled_sigma_chi.pract_dc;
					pract_dc_bin_pair_sync_bin_shuffled_sigma_chi(i,j) = practCE.bin_pair_sync_bin_shuffled_sigma_chi.pract_dc;
					
					pract_dc_bin_shuffled_phase_bin_sigma_chi(i,j) = practCE.bin_shuffled_phase_bin_sigma_chi.pract_dc;
					pract_dc_bin_shuffled_phase_bin_mean_pair_sync(i,j) = practCE.bin_shuffled_phase_bin_mean_pair_sync.pract_dc;
					
					pract_dc_bin_shuffled_phase_bin_shuffled_sigma_chi(i,j) = practCE.bin_phase_bin_shuffled_sigma_chi.pract_dc;
					
					% extract practical CD, and store it in arrays for different combinations of micro and macro variables
					pract_cd_bin_phase_bin_sigma_chi(i,j) = practCE.bin_phase_bin_sigma_chi.pract_cd;
					pract_cd_bin_phase_bin_mean_pair_sync(i,j) = practCE.bin_phase_bin_mean_pair_sync.pract_cd;
					pract_cd_bin_raw_signal_bin_sigma_chi(i,j) = practCE.bin_raw_signal_bin_sigma_chi.pract_cd;
					pract_cd_bin_raw_signal_bin_mean_pair_sync(i,j) = practCE.bin_raw_signal_bin_mean_pair_sync.pract_cd;
					pract_cd_bin_synchrony_bin_sigma_chi(i,j) = practCE.bin_synchrony_bin_sigma_chi.pract_cd;
					pract_cd_bin_synchrony_bin_mean_pair_sync(i,j) = practCE.bin_synchrony_bin_mean_pair_sync.pract_cd;
					pract_cd_bin_pair_sync_bin_sigma_chi(i,j) = practCE.bin_pair_sync_bin_sigma_chi.pract_cd;
					pract_cd_bin_pair_sync_bin_mean_pair_sync(i,j) = practCE.bin_pair_sync_bin_mean_pair_sync.pract_cd;
					
					% shuffled micro and macro variables
					pract_cd_bin_phase_bin_shuffled_sigma_chi(i,j) = practCE.bin_phase_bin_shuffled_sigma_chi.pract_cd;
					pract_cd_bin_raw_signal_bin_shuffled_sigma_chi(i,j) = practCE.bin_raw_signal_bin_shuffled_sigma_chi.pract_cd;
					pract_cd_bin_synchrony_bin_shuffled_sigma_chi(i,j) = practCE.bin_synchrony_bin_shuffled_sigma_chi.pract_cd;
					pract_cd_bin_pair_sync_bin_shuffled_sigma_chi(i,j) = practCE.bin_pair_sync_bin_shuffled_sigma_chi.pract_cd;
					
					pract_cd_bin_shuffled_phase_bin_sigma_chi(i,j) = practCE.bin_shuffled_phase_bin_sigma_chi.pract_cd;
					pract_cd_bin_shuffled_phase_bin_mean_pair_sync(i,j) = practCE.bin_shuffled_phase_bin_mean_pair_sync.pract_cd;
					
					pract_cd_bin_shuffled_phase_bin_shuffled_sigma_chi(i,j) = practCE.bin_phase_bin_shuffled_sigma_chi.pract_cd;
					
					clear practCE;
					clear macro_variables;
					clear micro_variables;
					
				end
			end
			
			% store practical CE, DC, & CE for all combinations of micro and macro variables and all parameter combinations 
			% of beta and A in one struct
			
			% practical CE
			practCE.pract_ce_bin_phase_bin_sigma_chi = pract_ce_bin_phase_bin_sigma_chi;
			practCE.pract_ce_bin_phase_bin_mean_pair_sync = pract_ce_bin_phase_bin_mean_pair_sync;
			practCE.pract_ce_bin_raw_signal_bin_sigma_chi = pract_ce_bin_raw_signal_bin_sigma_chi;
			practCE.pract_ce_bin_raw_signal_bin_mean_pair_sync = pract_ce_bin_raw_signal_bin_mean_pair_sync;
			practCE.pract_ce_bin_synchrony_bin_sigma_chi = pract_ce_bin_synchrony_bin_sigma_chi;
			practCE.pract_ce_bin_synchrony_bin_mean_pair_sync = pract_ce_bin_synchrony_bin_mean_pair_sync;
			practCE.pract_ce_bin_pair_sync_bin_sigma_chi = pract_ce_bin_pair_sync_bin_sigma_chi;
			practCE.pract_ce_bin_pair_sync_bin_mean_pair_sync = pract_ce_bin_pair_sync_bin_mean_pair_sync;
			
			% shuffled micro and macro variables
			practCE.pract_ce_bin_phase_bin_shuffled_sigma_chi = pract_ce_bin_phase_bin_shuffled_sigma_chi;
			practCE.pract_ce_bin_raw_signal_bin_shuffled_sigma_chi = pract_ce_bin_raw_signal_bin_shuffled_sigma_chi;
			practCE.pract_ce_bin_synchrony_bin_shuffled_sigma_chi = pract_ce_bin_synchrony_bin_shuffled_sigma_chi;
			practCE.pract_ce_bin_pair_sync_bin_shuffled_sigma_chi = pract_ce_bin_pair_sync_bin_shuffled_sigma_chi;
			
			practCE.pract_ce_bin_shuffled_phase_bin_sigma_chi = pract_ce_bin_shuffled_phase_bin_sigma_chi;
			practCE.pract_ce_bin_shuffled_phase_bin_mean_pair_sync = pract_ce_bin_shuffled_phase_bin_mean_pair_sync;
			
			practCE.pract_ce_bin_shuffled_phase_bin_shuffled_sigma_chi = pract_ce_bin_phase_bin_shuffled_sigma_chi;
			
			% practical DC
			practCE.pract_dc_bin_phase_bin_sigma_chi = pract_dc_bin_phase_bin_sigma_chi;
			practCE.pract_dc_bin_phase_bin_mean_pair_sync = pract_dc_bin_phase_bin_mean_pair_sync;
			practCE.pract_dc_bin_raw_signal_bin_sigma_chi = pract_dc_bin_raw_signal_bin_sigma_chi;
			practCE.pract_dc_bin_raw_signal_bin_mean_pair_sync = pract_dc_bin_raw_signal_bin_mean_pair_sync;
			practCE.pract_dc_bin_synchrony_bin_sigma_chi = pract_dc_bin_synchrony_bin_sigma_chi;
			practCE.pract_dc_bin_synchrony_bin_mean_pair_sync = pract_dc_bin_synchrony_bin_mean_pair_sync;
			practCE.pract_dc_bin_pair_sync_bin_sigma_chi = pract_dc_bin_pair_sync_bin_sigma_chi;
			practCE.pract_dc_bin_pair_sync_bin_mean_pair_sync = pract_dc_bin_pair_sync_bin_mean_pair_sync;
			
			% shuffled micro and macro variables
			practCE.pract_dc_bin_phase_bin_shuffled_sigma_chi = pract_dc_bin_phase_bin_shuffled_sigma_chi;
			practCE.pract_dc_bin_raw_signal_bin_shuffled_sigma_chi = pract_dc_bin_raw_signal_bin_shuffled_sigma_chi;
			practCE.pract_dc_bin_synchrony_bin_shuffled_sigma_chi = pract_dc_bin_synchrony_bin_shuffled_sigma_chi;
			practCE.pract_dc_bin_pair_sync_bin_shuffled_sigma_chi = pract_dc_bin_pair_sync_bin_shuffled_sigma_chi;
			
			practCE.pract_dc_bin_shuffled_phase_bin_sigma_chi = pract_dc_bin_shuffled_phase_bin_sigma_chi;
			practCE.pract_dc_bin_shuffled_phase_bin_mean_pair_sync = pract_dc_bin_shuffled_phase_bin_mean_pair_sync;
			
			practCE.pract_dc_bin_shuffled_phase_bin_shuffled_sigma_chi = pract_dc_bin_phase_bin_shuffled_sigma_chi;
			
			% practical CD
			practCE.pract_cd_bin_phase_bin_sigma_chi = pract_cd_bin_phase_bin_sigma_chi;
			practCE.pract_cd_bin_phase_bin_mean_pair_sync = pract_cd_bin_phase_bin_mean_pair_sync;
			practCE.pract_cd_bin_raw_signal_bin_sigma_chi = pract_cd_bin_raw_signal_bin_sigma_chi;
			practCE.pract_cd_bin_raw_signal_bin_mean_pair_sync = pract_cd_bin_raw_signal_bin_mean_pair_sync;
			practCE.pract_cd_bin_synchrony_bin_sigma_chi = pract_cd_bin_synchrony_bin_sigma_chi;
			practCE.pract_cd_bin_synchrony_bin_mean_pair_sync = pract_cd_bin_synchrony_bin_mean_pair_sync;
			practCE.pract_cd_bin_pair_sync_bin_sigma_chi = pract_cd_bin_pair_sync_bin_sigma_chi;
			practCE.pract_cd_bin_pair_sync_bin_mean_pair_sync = pract_cd_bin_pair_sync_bin_mean_pair_sync;
			
			% shuffled micro and macro variables
			practCE.pract_cd_bin_phase_bin_shuffled_sigma_chi = pract_cd_bin_phase_bin_shuffled_sigma_chi;
			practCE.pract_cd_bin_raw_signal_bin_shuffled_sigma_chi = pract_cd_bin_raw_signal_bin_shuffled_sigma_chi;
			practCE.pract_cd_bin_synchrony_bin_shuffled_sigma_chi = pract_cd_bin_synchrony_bin_shuffled_sigma_chi;
			practCE.pract_cd_bin_pair_sync_bin_shuffled_sigma_chi = pract_cd_bin_pair_sync_bin_shuffled_sigma_chi;
			
			practCE.pract_cd_bin_shuffled_phase_bin_sigma_chi = pract_cd_bin_shuffled_phase_bin_sigma_chi;
			practCE.pract_cd_bin_shuffled_phase_bin_mean_pair_sync = pract_cd_bin_shuffled_phase_bin_mean_pair_sync;
			
			practCE.pract_cd_bin_shuffled_phase_bin_shuffled_sigma_chi = pract_cd_bin_phase_bin_shuffled_sigma_chi;
			
			
			%% storing practical measures for different micro & macro variables
			
			% saved filenames consist of
			% network name + type of causal emergence + number of datapoints + time-lag
			save([pathout_data_pract_ce network '_practCE_' num2str(npoints) '_' num2str(tau) '.mat'], ...
				'practCE');
			
			clear pract_ce_bin_phase_bin_sigma_chi;
			clear pract_ce_bin_phase_bin_mean_pair_sync;
			clear pract_ce_bin_raw_signal_bin_sigma_chi;
			clear pract_ce_bin_raw_signal_bin_mean_pair_sync;
			clear pract_ce_bin_synchrony_bin_sigma_chi;
			clear pract_ce_bin_synchrony_bin_mean_pair_sync;
			clear pract_ce_bin_pair_sync_bin_sigma_chi;
			clear pract_ce_bin_pair_sync_bin_mean_pair_sync;
			
			clear pract_ce_bin_phase_bin_shuffled_sigma_chi;
			clear pract_ce_bin_raw_signal_bin_shuffled_sigma_chi;
			clear pract_ce_bin_synchrony_bin_shuffled_sigma_chi;
			clear pract_ce_bin_pair_sync_bin_shuffled_sigma_chi;
			clear pract_ce_bin_shuffled_phase_bin_sigma_chi;
			clear pract_ce_bin_shuffled_phase_bin_mean_pair_sync;
			clear pract_ce_bin_shuffled_phase_bin_shuffled_sigma_chi;
			
			clear pract_dc_bin_phase_bin_sigma_chi;
			clear pract_dc_bin_phase_bin_mean_pair_sync;
			clear pract_dc_bin_raw_signal_bin_sigma_chi;
			clear pract_dc_bin_raw_signal_bin_mean_pair_sync;
			clear pract_dc_bin_synchrony_bin_sigma_chi;
			clear pract_dc_bin_synchrony_bin_mean_pair_sync;
			clear pract_dc_bin_pair_sync_bin_sigma_chi;
			clear pract_dc_bin_pair_sync_bin_mean_pair_sync;
			
			clear pract_dc_bin_phase_bin_shuffled_sigma_chi;
			clear pract_dc_bin_raw_signal_bin_shuffled_sigma_chi;
			clear pract_dc_bin_synchrony_bin_shuffled_sigma_chi;
			clear pract_dc_bin_pair_sync_bin_shuffled_sigma_chi;
			clear pract_dc_bin_shuffled_phase_bin_sigma_chi;
			clear pract_dc_bin_shuffled_phase_bin_mean_pair_sync;
			clear pract_dc_bin_shuffled_phase_bin_shuffled_sigma_chi;
			
			clear pract_cd_bin_phase_bin_sigma_chi;
			clear pract_cd_bin_phase_bin_mean_pair_sync;
			clear pract_cd_bin_raw_signal_bin_sigma_chi;
			clear pract_cd_bin_raw_signal_bin_mean_pair_sync;
			clear pract_cd_bin_synchrony_bin_sigma_chi;
			clear pract_cd_bin_synchrony_bin_mean_pair_sync;
			clear pract_cd_bin_pair_sync_bin_sigma_chi;
			clear pract_cd_bin_pair_sync_bin_mean_pair_sync;
			
			clear pract_cd_bin_phase_bin_shuffled_sigma_chi;
			clear pract_cd_bin_raw_signal_bin_shuffled_sigma_chi;
			clear pract_cd_bin_synchrony_bin_shuffled_sigma_chi;
			clear pract_cd_bin_pair_sync_bin_shuffled_sigma_chi;
			clear pract_cd_bin_shuffled_phase_bin_sigma_chi;
			clear pract_cd_bin_shuffled_phase_bin_mean_pair_sync;
			clear pract_cd_bin_shuffled_phase_bin_shuffled_sigma_chi;
			
		end
		
		clear practCE;
		
	end

end 