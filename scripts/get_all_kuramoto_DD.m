function get_all_kuramoto_DD(network, A_vec, beta_vec, coupling_matrices, taus, tau_steps, all_npoints, ...
		method, pathout_data_sim_time_series, pathout_data_dd, kraskov_param)

	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		for z = 1:length(taus);
			tau = taus(z);
			
			rng(1);
			for i = 1:size(coupling_matrices, 3);
				A_str = param2str(A_vec(i));
				
				for j = 1:length(beta_vec)
					beta_str = param2str(beta_vec(j));
					
					% load simulated model with given A and beta
					load([pathout_data_sim_time_series network '_phase_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'phase');
					load([pathout_data_sim_time_series network '_raw_signal_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'raw_signal');
					load([pathout_data_sim_time_series network '_sigma_chi_' A_str '_' beta_str '_' num2str(npoints)  '.mat'], ...
						'sigma_chi');
					load([pathout_data_sim_time_series network '_synchrony_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'synchrony');
					load([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'mean_pair_sync');
					load([pathout_data_sim_time_series network '_pair_sync_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
						'pair_sync');
				
					% micro variables: binarized phase, raw signal, synchrony, pairwise synchrony
					% macro variables: binarized chimera-index (sigma_chi), average pairwise synchrony
					
					% ---------------------------------------------------------------------------------------------------------------------------------------
					% dynamical dependence
					% ---------------------------------------------------------------------------------------------------------------------------------------
					
					% store micro and macro variables in two structs (variables are transposed, so have time-points in rows)
					macro_variables.sigma_chi = sigma_chi';
					macro_variables.mean_pair_sync = mean_pair_sync';

					micro_variables.phase = phase';
					micro_variables.raw_signal = raw_signal';
					micro_variables.synchrony = synchrony';
					micro_variables.pair_sync = pair_sync';
					
					% calculate practical CE
					
					% get practical CE, DC, & CE for all combinations of micro and macro variables
					DD = get_DD(micro_variables, macro_variables, method, tau, tau_steps, kraskov_param);
					
					% extract DD for different combinations of micro and macro variables, and store it in arrays for 
					% different parameter combinations of beta and A
					dd_phase_sigma_chi(i,j) = DD.phase_sigma_chi;
					dd_phase_mean_pair_sync(i,j) = DD.phase_mean_pair_sync;
					dd_raw_signal_sigma_chi(i,j) = DD.raw_signal_sigma_chi;
					dd_raw_signal_mean_pair_sync(i,j) = DD.raw_signal_mean_pair_sync;
					dd_synchrony_sigma_chi(i,j) = DD.synchrony_sigma_chi;
					dd_synchrony_mean_pair_sync(i,j) = DD.synchrony_mean_pair_sync;
					dd_pair_sync_sigma_chi(i,j) = DD.pair_sync_sigma_chi;
					dd_pair_sync_mean_pair_sync(i,j) = DD.pair_sync_mean_pair_sync;
					
					clear DD;
					clear macro_variables;
					clear micro_variables;
					
					% store micro and macro variables in two structs (variables are transposed, so have time-points in rows)
					macro_variables.synchrony = synchrony';
					macro_variables.pair_sync = pair_sync';

					micro_variables.phase = phase';
					micro_variables.raw_signal = raw_signal';
					
					% calculate practical CE
					
					% get practical CE, DC, & CE for all combinations of micro and macro variables
					DD = get_DD(micro_variables, macro_variables, method, tau, tau_steps, kraskov_param);
					
					% extract DD for different combinations of micro and macro variables, and store it in arrays for 
					% different parameter combinations of beta and A
					dd_phase_synchrony(i,j) = DD.phase_synchrony;
					dd_phase_pair_sync(i,j) = DD.phase_pair_sync;
					dd_raw_signal_synchrony(i,j) = DD.raw_signal_synchrony;
					dd_raw_signal_pair_sync(i,j) = DD.raw_signal_pair_sync;
					
					clear DD;
					clear macro_variables;
					clear micro_variables;
					
				end
			end
			
			% store practical CE, DC, & CE for all combinations of micro and macro variables and all parameter combinations 
			% of beta and A in one struct
			
			% dynamical dependence - top-level macro variables
			DD.dd_phase_sigma_chi = dd_phase_sigma_chi;
			DD.dd_phase_mean_pair_sync = dd_phase_mean_pair_sync;
			DD.dd_raw_signal_sigma_chi = dd_raw_signal_sigma_chi;
			DD.dd_raw_signal_mean_pair_sync = dd_raw_signal_mean_pair_sync;
			DD.dd_synchrony_sigma_chi = dd_synchrony_sigma_chi;
			DD.dd_synchrony_mean_pair_sync = dd_synchrony_mean_pair_sync;
			DD.dd_pair_sync_sigma_chi = dd_pair_sync_sigma_chi;
			DD.dd_pair_sync_mean_pair_sync = dd_pair_sync_mean_pair_sync;
			
			% mid-level macro variables
			DD.dd_phase_synchrony = dd_phase_synchrony;
			DD.dd_phase_pair_sync = dd_phase_pair_sync;
			DD.dd_raw_signal_synchrony = dd_raw_signal_synchrony;
			DD.dd_raw_signal_pair_sync = dd_raw_signal_pair_sync;
			
			%% storing practical measures for different micro & macro variables
			
			% saved filenames consist of
			% network name + type of causal emergence + number of datapoints + time-lag
			save([pathout_data_dd network '_DD_' num2str(npoints) '_' num2str(tau) '.mat'], ...
				'DD');
			
			clear dd_phase_sigma_chi;
			clear dd_phase_mean_pair_sync;
			clear dd_raw_signal_sigma_chi;
			clear dd_raw_signal_mean_pair_sync;
			clear dd_synchrony_sigma_chi;
			clear dd_synchrony_mean_pair_sync;
			clear dd_pair_sync_sigma_chi;
			clear dd_pair_sync_mean_pair_sync;
			clear dd_phase_synchrony;
			clear dd_phase_pair_sync;
			clear dd_raw_signal_synchrony;
			clear dd_raw_signal_pair_sync;
			
		end
		
		clear DD;
		
	end
