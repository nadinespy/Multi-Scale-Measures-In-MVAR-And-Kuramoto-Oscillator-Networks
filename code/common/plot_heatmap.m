% heatmaps using imagesc()

function plot_heatmap(data, file_names, titles, x_axis, y_axis, x_label, y_label, pathout_plots) 
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'data', @iscell);
	addRequired(p,'file_names', @iscell)
	addRequired(p,'titles', @iscell);
	addRequired(p,'x_axis', @iscell);
	addRequired(p,'y_axis', @iscell);
	addRequired(p,'x_label', @ischar);
	addRequired(p,'y_label', @ischar);
	addRequired(p,'pathout_plots', @ischar);

	parse(p, data, file_names, titles, x_axis, y_axis, x_label, ...
		y_label, pathout_plots);
	
	data			= p.Results.data;
	file_names		= p.Results.file_names;
	titles			= p.Results.titles;
	x_axis			= p.Results.x_axis;
	y_axis			= p.Results.y_axis;
	x_label			= p.Results.x_label;
	y_label			= p.Results.y_label;
	pathout_plots	= p.Results.pathout_plots;
	
	for i = 1:length(data)
		
		convert_data = data{i};
		
		if istable(data{i})
			convert_data = table2array(convert_data);
		end
		
		fig = figure;
		
		imagesc(convert_data);
		colormap(bluewhitered);
		colorbar;
		
		hColorbar = colorbar;
		
		xticks(1:(size(x_axis, 2)));
		yticks(1:(size(y_axis, 2)));
		
		set(gca,'TickLength',[0 0])
		yticklabels(y_axis);
		xticklabels(x_axis);
		
		ylabel(y_label);
		xlabel(x_label);
		
		title(titles{i});
		
		location = string(strcat(pathout_plots, file_names{i}, '.png'));
		exportgraphics(gcf, location);
		
		close all;
	end
	
end 