function get_km_met_chi_heatmaps(network, npoints, x_axis, y_axis, x_label, y_label, ...
		pathout_data_sync, pathout_plots_chi, pathout_plots_met)

	set(0,'DefaultFigureVisible','off');
	
	for q = 1:length(npoints);
		
		load([pathout_data_sync network '_mean_chi_' num2str(npoints(q)) '.mat'], ...
				'mean_chi');
		load([pathout_data_sync network '_mean_met_' num2str(npoints(q)) '.mat'], ...
				'mean_met');
			
		variables = {mean_chi};	
		file_names = {'mean_chi'};
		titles = {{['mean chi, ' 'npoints = ' num2str(npoints(q))]}};
			
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints(q), ...
				pathout_plots_chi);
			
		variables = {mean_met};	
		file_names = {'mean_met'};
		titles = {{['mean met, ' 'npoints = ' num2str(npoints(q))]}};
			
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints(q), ...
				pathout_plots_met);
			
	end 
	
	set(0,'DefaultFigureVisible','off');
				
end 
