function get_all_kuramoto_met_chi(network, A_vec, all_npoints, ...
		pathout_data_sync)
	
	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		for p = 1:length(A_vec);
			A_str = param2str(A_vec(p));
			
			% load synchronies for a given value of A and across beta;
			% rows: betas; columns: communities; 3rd dimension: time-points

			load([pathout_data_sync network '_synchronies_'  A_str '_' num2str(npoints) '.mat'], ...
				'synchronies');
			
			for k = 1:size(synchronies, 1);
				
				% calculate metastability:
				% mean of variance over time for each synchrony
				sigma_met = squeeze(synchronies(k,:,:));
				sigma_met_mean(p,k) = mean(var(sigma_met')); % sigma_met': time x synchronies
				
				
				% calculate chimera states
				% mean of variance over synchronies for each time-point
				sigma_chi = squeeze(synchronies(k,:,:));
				sigma_chi_mean(p,k) = mean(var(sigma_chi));  % sigma_chi: synchronies x time

			end		
			
		end 
		
		save([pathout_data_sync network '_sigma_chi_mean_' num2str(npoints) '.mat'], ...
				'sigma_chi_mean');
		save([pathout_data_sync network '_sigma_met_mean_' num2str(npoints) '.mat'], ...
				'sigma_met_mean');
	end 
end 

			
			