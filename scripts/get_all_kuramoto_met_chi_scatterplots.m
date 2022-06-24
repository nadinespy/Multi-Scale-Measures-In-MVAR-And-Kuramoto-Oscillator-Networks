function get_all_kuramoto_met_chi_scatterplots(network, A_vec, beta_vec, all_npoints, x_label, y_label, ...
		pathout_data_sync, pathout_plots_sigma_chi, pathout_plots_sigma_met)

	set(0,'DefaultFigureVisible','off');
	
	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		n_A = length(A_vec);
		A = [1:3:n_A];
		
		for p = 1:length(A);
	
			g = A(p);
			A_str = param2str(A_vec(g));
			
			% load synchronies for a given value of A
			load([pathout_data_sync network '_synchronies_'  A_str '_' num2str(npoints) '.mat'], ...
				'synchronies');
			
			% calculate metastability
			sigma_met_mean = [];
			for i = 1:length(beta_vec);
				sigma_met = squeeze(synchronies(i,:,:));
				sigma_met_mean(i) = mean(var(sigma_met')); % sigma_met': time x synchronies
			end							       % mean of variance over time for each synchrony
			
			% calculate chimera states
			sigma_chi_mean = [];
			for i = 1:length(beta_vec);
				sigma_chi = squeeze(synchronies(i,:,:));
				sigma_chi_mean(i) = mean(var(sigma_chi));  % sigma_chi: synchronies x time
			end								 % mean of variance over synchronies for each time-point
			
			
			variables = {sigma_chi_mean, sigma_met_mean};
			
			file_names = {'sigma_chi_mean', ...
				'sigma_met_mean'};
			
			titles = {{['sigma chi mean' ', npoints = ' num2str(npoints)]}, ...
				{['sigma met mean' ', npoints = ' num2str(npoints)]}};
			
			pathout_plots = {pathout_plots_sigma_chi, pathout_plots_sigma_met};
			
			plot_scatterplot_met_chi(variables, A_vec(g), beta_vec, file_names, titles, x_label, y_label, ...
				network, npoints, pathout_plots) 
		end 
	end 
end 
