% heatmaps using imagesc()

function plot_heatmap(data, file_names, titles, x_axis, y_axis, x_label, y_label, pathout_plots) 
	
	for i = 1:length(data)
			
		convert_data = table2array(data{i});
		fig = figure;
		
		imagesc(convert_data);
		colormap(bluewhitered);
		colorbar;
		
		hColorbar = colorbar;
		
		xticks(1:size(x_axis, 2));
		yticks(1:size(y_axis, 2));
		
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