function get_km_met_chi_heatmaps(network, metastability_results, time_lengths, ...
		x_axis, y_axis, x_label, y_label, pathout_data_metastability, ...
		pathout_plots_metastability)
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'metastability_results', @iscell)
	addRequired(p,'time_lengths', @isdouble);
	addRequired(p,'x_axis', @iscell);
	addRequired(p,'y_axis', @iscell);
	addRequired(p,'x_label', @ischar);
	addRequired(p,'y_label', @ischar);
	addRequired(p,'pathout_data_metastability', @ischar);
	addRequired(p,'pathout_plots_metastability', @ischar);

	parse(p, network, metastability_results, time_lengths, x_axis, ...
		y_axis, x_label, y_label, pathout_data_metastability, ...
		pathout_plots_metastability);
	
	network				= p.Results.network;
	metastability_results		= p.Results.metastability_results;
	time_lengths			= p.Results.time_lengths;
	x_axis				= p.Results.x_axis;
	y_axis				= p.Results.y_axis;
	x_label				= p.Results.x_label;
	y_label				= p.Results.y_label;
	pathout_data_metastability	= p.Results.pathout_data_metastability;
	pathout_plots_metastability	= p.Results.pathout_plots_metastability;

	set(0,'DefaultFigureVisible','off');
	
	for q = 1:length(time_lengths);
		
		for k = 1:length(metastability_results)
			
			meta_var = load([pathout_data_metastability network '_' metastability_results{k} '_' ...
				num2str(time_lengths(q)) '.mat'], metastability_results{k});
			
			meta_var_fieldnames = fieldnames(meta_var);
			meta_var = {meta_var.(meta_var_fieldnames{1})};
			
			% flip rows in order
			meta_var{1} = flip(meta_var{1},1);
			
			file_names = {{[network '_' meta_var_fieldnames{1} '_' num2str(time_lengths(q))]}};
			titles = {[meta_var_fieldnames{1} ', time_lengths = ' num2str(time_lengths(q))]};
			
			set(0,'DefaultFigureVisible','on');
			
			plot_heatmap(meta_var, file_names, titles, x_axis, y_axis, x_label, y_label, ...
				pathout_plots_metastability);
			
		end
			
	end 
				
end 
