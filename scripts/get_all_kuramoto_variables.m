function get_all_kuramoto_variables(network, intra_comm_size, n_communities, coupling_matrices, ...
	param1_vec, param2_vec, param3_vec, pathout_data_sim_time_series, pathout_data_sync);

	% sim_kuramoto_oscillators() obtains variables 'phase', 'sigma_chi', and 'synchrony';
	% synchronies for different values of beta and given value of A are stored and saved in 'synchronies';
	% 'grand_mean_pair_sync' and 'raw_values' are derived using 'synchrony' and 'phase', respectively

	% SIMULATE KURAMOTO OSCILLATORS: outputs of sim_kuramoto_oscillators() are
	%	- phases				(MICRO)
	%	- chimera-index			(MACRO)
	%	- synchronies			(MICRO)

	% GET FURTHER MICRO AND MACRO VARIABLES:
	%	- raw signal (cos(phase))	(MICRO)
	%	- average pairwise synchrony	(MACRO)

	for q = 1:length(param1_vec);
		param1_str = num2str(param1_vec(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			param2_str = param2str(param2_vec(i));
		
			for j = 1:length(param3_vec)
				param3_str = param2str(param3_vec);
				
				fprintf('get_all_kuramoto_variables - loop indices: npoints: %d, A: %d, beta: %d\n', q, i, j);
			
				[phase, sigma_chi, synchrony] = sim_kuramoto_oscillators(intra_comm_size, n_communities, ...
					coupling_matrix, param3_vec(j), param1_vec(q));
			
				% store synchronies for a given A, and across beta;
				% rows: betas; columns: communities; 3rd dimension: time-points
				synchronies(j,:,:) = synchrony;
			
				% MICRO VARIABLES:
				%	- phases,
				%	- raw signal,
				%	- synchronies
				%	- pair_sync
				%	- shuffled phases as a control
				
				% raw signal
				raw_signal = cos(phase);

				% MACRO VARIABLES for practical measures for causal emergence:
				%	- synchronies
				%	- pairwise synchrony
				%	- variance of synchronies (sigma_chi)
				%	- mean pairwise synchrony between communities (mean_pair_sync)
				%	- shuffled sigma_chi as a control
				
				% global mean pairwise synchrony
				[pair_sync, mean_pair_sync] = get_kuramoto_pair_sync(n_communities, synchrony, param1_vec(q));

				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				save([pathout_data_sim_time_series network '_phase_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'phase');
				save([pathout_data_sim_time_series network '_sigma_chi_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'sigma_chi');
				save([pathout_data_sim_time_series network '_synchrony_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'synchrony');
				save([pathout_data_sim_time_series network '_pair_sync_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'pair_sync');
				save([pathout_data_sim_time_series network '_mean_pair_sync_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'mean_pair_sync');
				save([pathout_data_sim_time_series network '_raw_signal_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'raw_signal');
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout_data_sync network '_synchronies_'  param2_str '_' param1_str '.mat'], ...
				'synchronies');
			
			clear synchrony;
			clear synchronies;
		end
	end
end 

