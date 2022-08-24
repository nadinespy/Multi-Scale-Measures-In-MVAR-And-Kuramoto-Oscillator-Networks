function get_km_met_chi(network, npoints, model_params, ...
		pathout_data_sync)
	
	model_param1 = model_params{1};
	
	for q = 1:length(npoints);
		
		for p = 1:length(model_param1);
			model_param1_str = param2str(model_param1(p));
			
			% load synchronies for a given value of model_param1 (A_vec) and across model_params{2} (beta_vec);
			% rows: betas; columns: communities; 3rd dimension: time-points

			load([pathout_data_sync network '_synchronies_'  model_param1_str '_' num2str(npoints(q)) '.mat'], ...
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
		
		save([pathout_data_sync network '_mean_chi_' num2str(npoints(q)) '.mat'], ...
				'mean_chi');
		save([pathout_data_sync network '_mean_met_' num2str(npoints(q)) '.mat'], ...
				'mean_met');
	end 
end 

			
			
