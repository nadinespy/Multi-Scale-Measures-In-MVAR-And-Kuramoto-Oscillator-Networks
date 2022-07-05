function get_all_practCE_heatmaps(network, time_series_length, measure_params, x_axis, y_axis, x_label, y_label, ...
		pathout_data_pract_ce, pathout_plots_pract_ce)

	set(0,'DefaultFigureVisible','off');
	
	measure_param1 = measure_params{1};
	
	for q = 1:length(time_series_length);
		
		for z = 1:length(measure_param1);
			
			load([pathout_data_pract_ce network '_all_practCE_' num2str(time_series_length(q)) '_' num2str(measure_param1(z)) '.mat'], ...
				'all_practCE');
			
			fieldnames_all_practCE = fieldnames(all_practCE);

			for k = 1:length(fieldnames(all_practCE));

				variables{k} = all_practCE.(fieldnames_all_practCE{k,1});
				% file names consist of network + practical measure + micro variable + macro variable + ...
				% disc number/pract CE method + time-series length + measure param1
				file_names{k} = fieldnames_all_practCE{k,1};
				titles{k} = {[fieldnames_all_practCE{k,1}] ['npoints = ' num2str(time_series_length(q)) ', tau = ' num2str(measure_param1(z))]};
				
			end
		
			plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, time_series_length(q), ...
				pathout_plots_pract_ce, measure_param1(z));
		
		end
	end
	
end 
