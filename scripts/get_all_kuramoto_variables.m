function get_all_kuramoto_variables(network, model_params, coupling_matrices, time_series_length, ...
		pathout_data_sim_time_series, pathout_data_sync);

	model_param1 = model_params{1};
	model_param2 = model_params{2};
	model_param3 = model_params{3};
	model_param4 = model_params{4};
	
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

	for q = 1:length(time_series_length);
		time_series_length_str = num2str(time_series_length(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			model_param1_str = param2str(model_param1(i));
		
			for j = 1:length(model_param2)
				model_param2_str = param2str(model_param2(j));
				
				fprintf('get_all_kuramoto_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d\n', q, i, j);
			
				[phase, chi, sync] = sim_kuramoto_oscillators(model_param3, model_param4, ...
					coupling_matrix, model_param2(j), time_series_length(q));
			
				% store synchronies for a given A, and across beta;
				% rows: betas; columns: communities; 3rd dimension: time-points
				synchronies(j,:,:) = sync;
			
				% MICRO VARIABLES:
				%	- phases,
				%	- raw signal,
				%	- synchronies
				%	- pair_sync
				%	- shuffled phases as a control
				
				% raw signal
				raw = cos(phase);

				% MACRO VARIABLES for practical measures for causal emergence:
				%	- synchronies
				%	- pairwise synchrony
				%	- variance of synchronies (sigma_chi)
				%	- mean pairwise synchrony between communities (mean_pair_sync)
				%	- shuffled sigma_chi as a control
				
				% global mean pairwise synchrony
				[p_sync, mp_sync] = get_kuramoto_pair_sync(model_param4, sync, time_series_length(q));

				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				save([pathout_data_sim_time_series network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'phase');
				save([pathout_data_sim_time_series network '_chi_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'chi');
				save([pathout_data_sim_time_series network '_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'sync');
				save([pathout_data_sim_time_series network '_p_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'p_sync');
				save([pathout_data_sim_time_series network '_mp_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'mp_sync');
				save([pathout_data_sim_time_series network '_raw_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'raw');
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout_data_sync network '_synchronies_'  model_param1_str '_' time_series_length_str '.mat'], ...
				'synchronies');
			
			clear synchrony;
			clear synchronies;
		end
	end
end 

