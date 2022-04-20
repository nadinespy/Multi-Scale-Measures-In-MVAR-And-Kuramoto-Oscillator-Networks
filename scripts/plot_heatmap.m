% heatmaps using imagesc()

function plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots, tau) 
	
	for i = 1:size(variables,2)
			
			figure;
			
			imagesc(variables{i});
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
			
			if nargin < 11;
				exportgraphics(gcf, [pathout_plots network '_' file_names{i} '_' num2str(npoints) '.png']);
			else 
				exportgraphics(gcf, [pathout_plots network '_' file_names{i} '_' num2str(npoints) '_' num2str(tau) '.png']);
			end
			
			close all;
	end
	
end 