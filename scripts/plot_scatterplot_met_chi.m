function plot_scatterplot_met_chi(variables, A, beta_vec, file_names, titles, ...
		x_label, y_label, network, npoints, pathout_plots) 
	
	% coupling_vec(20) = 0.2131
	% A = [2, 4, 6, 8, 10];
	A_str = param2str(A);
	
	for i = 1:size(variables,2)

		variable_temp = variables{i};
		title_temp = [titles{i}, {['A = ' A]}];
		
		figure;
		scatter(beta_vec, variable_temp, 60, 'filled');
		title(title_temp)
		ylabel(y_label);
		xlabel(x_label);

		exportgraphics(gcf, [pathout_plots{i} network '_' file_names{i} '_' A_str '_' num2str(npoints) '.png']);
		
		close all;
		
	end
end 