function get_all_kuramoto_variables(network, intra_comm_size, n_communities, kuramoto_coupling_matrices, ...
	A_vec, beta_vec, all_npoints, pathout_data_sim_time_series, pathout_data_sync);

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

	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
	
		for i = 1:size(kuramoto_coupling_matrices, 3);
			kuramoto_coupling_matrix = kuramoto_coupling_matrices(:,:,i);
			A_str = param2str(A_vec(i));
		
			for j = 1:length(beta_vec)
				beta = beta_vec(j);
				beta_str = param2str(beta);
			
				[phase, sigma_chi, synchrony] = sim_kuramoto_oscillators(intra_comm_size, n_communities, ...
					kuramoto_coupling_matrix, beta, npoints);
			
				% store synchronies for a given A, and across beta;
				% rows: betas; columns: communities; 3rd dimension: time-points
				synchronies(j,:,:) = synchrony;
			
				% MICRO VARIABLES:
				%	- phases,
				%	- raw signal,
				%	- synchronies
				
				% raw signal
				raw_signal = cos(phase);
			
				% MACRO VARIABLES for practical measures for causal emergence:
				%	- variance of synchronies (sigma_chi)
				%	- global mean pairwise synchrony between communities (grand_mean_pair_sync)
				
				% global mean pairwise synchrony
				grand_mean_pair_sync = get_kuramoto_grand_mean_pair_sync(n_communities, synchrony, npoints);
				
				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				save([pathout_data_sim_time_series network '_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'phase');
				save([pathout_data_sim_time_series network '_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'sigma_chi');
				save([pathout_data_sim_time_series network '_synchrony_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'synchrony');
				save([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'grand_mean_pair_sync');
				save([pathout_data_sim_time_series network '_raw_signal_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'raw_signal');
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout_data_sync network '_synchronies_'  A_str '_' num2str(npoints) '.mat'], ...
				'synchronies');
			
			clear synchrony;
			clear synchronies;
		end
	end
end 

