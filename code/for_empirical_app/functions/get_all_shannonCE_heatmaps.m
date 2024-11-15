function get_all_heatmaps(emergence_struct, x_axis, y_axis, x_label, y_label, ...
		pathout_data_pract_ce, pathout_plots_pract_ce)

	set(0,'DefaultFigureVisible','off');
	
	for q = 1:length(emergence_struct);
	
		micro_macro_fieldnames = fieldnames(emergence_struct(1,q).results);
		
		for k = 1:length(micro_macro_fieldnames);
			
			time_length_str = num2str(emergence_struct(1,q).time_length);
			
			titles{k} = {[micro_macro_fieldnames{k,1}, ', npoints = ' num2str(emergence_struct(1,q).time_length)] ...
				[', measure = ', emergence_struct(1,q).measure, ...
				', method = ', emergence_struct(1,q).method, ...
				', timelag = ' num2str(emergence_struct(1,q).time_lag)]};
			
			if strcmp(lower(emergence_struct(1,q).method), 'discrete')
				titles_discrete{k} = {[', disc method = ' emergence_struct(1,q).disc_method, ...
					', bins = ' num2str(emergence_struct(1,q).bin_number)]};
				
				titles{1,k}{length(titles{1,k})+1} = titles_discrete{k};
			end
			
			if strcmp(lower(emergence_struct(1,q).method), 'kraskov')
				titles_kraskov{k} = {[', kraskov param = ' emergence_struct(1,q).kraskov_param]};
				
				titles{1,k}{length(titles{1,k})+1} = titles_kraskov{k};
			end
			
			if strcmp(lower(emergence_struct(1,q).measure), 'dd')
				titles_dd{k} = {[', time step = ' emergence_struct(1,q).time_step]};
				
				titles{1,k}{length(titles{1,k})+1} = titles_dd{k};
			end
			
			if strfind(emergence_struct(1,q).measure, 'phiid') == 1
				titles_phiid{k} = {[', red func = ' emergence_struct(1,q).red_func]};
				
				titles{1,k}{length(titles{1,k})+1} = titles_phiid{k};
			end
			
			file_names{k} = {[network '_' micro_macro_fieldnames{k,1} '_' emergence_struct(1,q).measure '_' ...
				lower(emergence_struct(1,q).method(1:4)) '_lag' num2str(emergence_struct(1,q).time_lag) ...
				'_length' time_length_str(1:end-3) 'k']};
			
		end
		
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, ...
			pathout_plots_pract_ce);
	end
end 
