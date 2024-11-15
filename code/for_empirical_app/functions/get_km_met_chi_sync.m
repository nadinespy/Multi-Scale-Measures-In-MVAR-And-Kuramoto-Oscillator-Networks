function get_km_met_chi_sync(network, model_sim_params, time_lengths, pathin, pathout)
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'model_sim_params', @isstruct);
	addRequired(p,'time_lengths', @isdouble);
	addRequired(p,'pathin', @isstruct);
	addRequired(p,'pathout', @isstruct);
	
	parse(p, network, model_sim_params, time_lengths, pathin, pathout);

	network		= p.Results.network;
	model_sim_params	= p.Results.model_sim_params;
	time_lengths	= p.Results.time_lengths;
	pathin		= p.Results.pathin;
	pathout		= p.Results.pathout;

	% extract cell arrays from structs
	model_sim_params_fieldnames	= fieldnames(model_sim_params);
	model_params1			= model_sim_params.(model_sim_params_fieldnames{1});
	model_params2			= model_sim_params.(model_sim_params_fieldnames{2});
	
	pathin_fieldnames			= fieldnames(pathin);
	pathout_data_sync			= pathin.(pathin_fieldnames{2});
	pathout_data_sim_time_series	= pathin.(pathin_fieldnames{1});
	
	pathout_fieldnames		= fieldnames(pathout);
	pathout_data_metastability	= pathout.(pathout_fieldnames{2});
	
	for q = 1:length(time_lengths);
		
		for p = 1:length(model_params1);
			model_param1_str = param2str(model_params1(p));
			
			% load synchronies for a given value of model_param1 (A_vec) and across model_params{2} (beta_vec);
			% rows: betas; columns: communities; 3rd dimension: time-points

			load([pathout_data_sync network '_synchronies_'  model_param1_str '_' num2str(time_lengths(q)) '.mat'], ...
				'synchronies');
			
			for k = 1:size(synchronies, 1);
				
				% calculate metastability:
				% mean of variance over time for each synchrony
				met = squeeze(synchronies(k,:,:));
				mean_met(p,k) = mean(var(met')); % mean_met': time x synchronies
				
				
				% calculate chimera states
				% mean of variance over synchronies for each time-point
				chi = squeeze(synchronies(k,:,:));
				mean_chi(p,k) = mean(var(chi));  % sigma_chi: synchronies x time

			end	
						
		end 
		
		save([pathout_data_metastability network '_mean_chi_' num2str(time_lengths(q)) '.mat'], ...
				'mean_chi');
		save([pathout_data_metastability network '_mean_met_' num2str(time_lengths(q)) '.mat'], ...
				'mean_met');
	end 
	
	for q = 1:length(time_lengths);
		
		for p = 1:length(model_params1);
			model_param1_str = param2str(model_params1(p));
			
			for j = 1:length(model_params2);
				model_param2_str = param2str(model_params2(j));
			
				load([pathout_data_sim_time_series network '_full_system_sync_'  model_param1_str '_' model_param2_str '_' num2str(time_lengths(q)) '.mat'], ...
					'full_system_sync');
			
				mean_full_system_sync(p,j) = mean(full_system_sync);
			end
		end
		
		save([pathout_data_metastability network '_mean_full_system_sync_' num2str(time_lengths(q)) '.mat'], ...
				'mean_full_system_sync');
	end
		
end 

			
			
