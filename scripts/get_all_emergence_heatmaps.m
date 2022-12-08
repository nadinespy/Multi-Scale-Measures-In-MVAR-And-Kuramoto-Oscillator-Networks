function get_all_emergence_heatmaps(network, emergence_struct, x_axis, y_axis, x_label, y_label, pathout_plots)

	set(0,'DefaultFigureVisible','off');
	
	o = 1;
	for q = 1:length(emergence_struct);
			
		micro_macro_fieldnames = fieldnames(emergence_struct(1,q).results);	
		micro_macro_fieldnames = fieldnames(emergence_struct(1,q).results);
		
		
		for k = 1:length(micro_macro_fieldnames);
			
			time_length_str = num2str(emergence_struct(1,q).time_length);
			data{o} = emergence_struct(1,q).results.(micro_macro_fieldnames{k});
			
			titles{1,o} = {[micro_macro_fieldnames{k,1}, ', npoints = ' num2str(emergence_struct(1,q).time_length), ', '] ...
				['measure = ', emergence_struct(1,q).measure, ', ', ...
				'method = ', emergence_struct(1,q).method, ', ', ...
				'timelag = ' num2str(emergence_struct(1,q).time_lag)]};
			
			if strcmp(lower(emergence_struct(1,q).method), 'discrete')
				titles_discrete = {['disc method = ' emergence_struct(1,q).disc_method, ', ', ...
					'bins = ' num2str(emergence_struct(1,q).bin_number)]};
				
				titles{1,o}(end+1) = titles_discrete;
			end
			
			if strcmp(lower(emergence_struct(1,q).method), 'kraskov')
				titles_kraskov = {['kraskov param = ' num2str(emergence_struct(1,q).kraskov_param)]};
				
				titles{1,o}(end+1) = titles_kraskov;
				
			end
			
			if strcmp(lower(emergence_struct(1,q).measure), 'dd')
				titles_dd = {['time step = ' num2str(emergence_struct(1,q).time_step)]};
				
				titles{1,o}(end+1) = titles_dd;
			end
			
			if strfind(emergence_struct(1,q).measure, 'phiid') == 1
				titles_phiid = {['red func = ' emergence_struct(1,q).red_func]};
				
				titles{1,o}(end+1) = titles_phiid;
			end
			
			file_names{o} = {[num2str(o) '_' network '_' micro_macro_fieldnames{k,1} '_' emergence_struct(1,q).measure '_' ...
				lower(emergence_struct(1,q).method(1:4)) '_lag' num2str(emergence_struct(1,q).time_lag) ...
				'_length' time_length_str(1:end-3) 'k']};
			
			titles{1,o}{1} = strrep(titles{1,o}{1},'_',' ');
			
			o = o+1;
		end
	end
	
	plot_heatmap(data, file_names, titles, x_axis, y_axis, x_label, y_label, ...
		pathout_plots);
end 
