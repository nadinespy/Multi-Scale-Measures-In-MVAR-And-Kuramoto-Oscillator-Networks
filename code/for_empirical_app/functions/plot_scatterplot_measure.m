function plot_scatterplot_measure(variables, A_vec, beta_vec, file_names, titles, x_label, y_label, network, ...
		npoints, pathout_plots, tau) 
	
	% coupling_vec(20) = 0.2131
	n_A = length(A_vec);
	A = [1:3:n_A];
	
	for i = 1:size(variables,2)
		
		for p = 1:length(A);
			
			g = A(p);
			A_str = param2str(A_vec(g));
			
			variable_temp = variables{i};
			title_temp = [titles{i}, {['A = ' num2str(A_vec(g))]}];
			
			figure;
			scatter(beta_vec, variable_temp(g,:), 60, 'filled');
			title(title_temp) 
			ylabel(y_label);
			xlabel(x_label);
			
			if nargin < 10;
				exportgraphics(gcf, [pathout_plots network '_' file_names{i} '_' A_str '_' num2str(npoints) '.png']);
			else
				exportgraphics(gcf, [pathout_plots network '_' file_names{i} '_' A_str '_' num2str(npoints) '_' num2str(tau) '.png']);
			end
			
			close all;
			
		end
		
	end
end 