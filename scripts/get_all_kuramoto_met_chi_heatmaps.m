function get_all_kuramoto_met_chi_heatmaps(network, all_npoints, x_axis, y_axis, x_label, y_label, ...
		pathout_data_sync, pathout_plots_sigma_chi, pathout_plots_sigma_met)

	set(0,'DefaultFigureVisible','off');
	
	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		load([pathout_data_sync network '_sigma_chi_mean_' num2str(npoints) '.mat'], ...
				'sigma_chi_mean');
		load([pathout_data_sync network '_sigma_met_mean_' num2str(npoints) '.mat'], ...
				'sigma_met_mean');
			
		variables = {sigma_chi_mean};	
		file_names = {'sigma_chi_mean'};
		titles = {{['sigma chi mean, ' 'npoints = ' num2str(npoints)]}};
			
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, ...
				pathout_plots_sigma_chi);
			
		variables = {sigma_met_mean};	
		file_names = {'sigma_met_mean'};
		titles = {{['sigma met mean, ' 'npoints = ' num2str(npoints)]}};
			
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, ...
				pathout_plots_sigma_met);
			
	end 
end 
