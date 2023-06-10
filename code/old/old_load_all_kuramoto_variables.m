function [phase, sigma_chi, synchrony, pair_sync, mean_pair_sync, raw_signal, shuffled_sigma_chi, shuffled_phase] = load_all_kuramoto_variables(network, A_vec, beta_vec, all_npoints, pathout_data_sim_time_series, pathout_data_sync);

	% load all variables generated in get_all_kuramoto_variables()

	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
	
		for i = 1:length(A_vec);
			A_str = param2str(A_vec(i));
		
			for j = 1:length(beta_vec)
				beta = beta_vec(j);
				beta_str = param2str(beta);
			
				% MICRO VARIABLES:
				%	- phases,
				%	- raw signal,
				%	- synchronies
				%	- pair_sync
				%	- shuffled phases as a control

				% MACRO VARIABLES for practical measures for causal emergence:
				%	- synchronies
				%	- pairwise synchrony
				%	- variance of synchronies (sigma_chi)
				%	- mean pairwise synchrony between communities (mean_pair_sync)
				%	- shuffled sigma_chi as a control

				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				phase = load([pathout_data_sim_time_series network '_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'phase');
				sigma_chi = load([pathout_data_sim_time_series network '_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'sigma_chi');
				synchrony = load([pathout_data_sim_time_series network '_synchrony_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'synchrony');
				pair_sync = load([pathout_data_sim_time_series network '_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'pair_sync');
				mean_pair_sync = load([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'mean_pair_sync');
				raw_signal = load([pathout_data_sim_time_series network '_raw_signal_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'raw_signal');
				shuffled_sigma_chi = load([pathout_data_sim_time_series network '_shuffled_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'shuffled_sigma_chi');
				shuffled_phase = load([pathout_data_sim_time_series network '_shuffled_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'shuffled_phase');
			end

		end
	end
end 

