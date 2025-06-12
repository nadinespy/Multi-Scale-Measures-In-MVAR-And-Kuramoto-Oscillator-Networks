% Function to create dual-axis plot
function plot_dual_lines(all_phase_lag, chi_met_data, dd_data, config, actual_coup_ratio, all_rmi, n_rmi_to_plot, network, pathout_plots_measures)
	figure;
	colors = lines(n_rmi_to_plot * 2);
	color_idx = 1;

	% Left y-axis (Chimera/Metastability)
	yyaxis left;
	for h = 1:n_rmi_to_plot
		plot(all_phase_lag, chi_met_data{config.data_idx}(h, :), 'LineWidth', 1.5, ...
			'DisplayName', config.name, 'LineStyle', '-', 'Color', colors(color_idx, :));
		hold on;
		color_idx = color_idx + 1;
	end
	ylabel(config.name, 'FontSize', 15, 'interpreter', 'latex');
	set(gca, 'YColor', 'k');

	% Add margin to left y-axis
	ylim_left = ylim;
	margin_left = 0.1 * (ylim_left(2) - ylim_left(1));
	ylim([ylim_left(1) - margin_left, ylim_left(2) + margin_left]);

	% Right y-axis (DD)
	yyaxis right;
	for h = 1:n_rmi_to_plot
		plot(all_phase_lag, dd_data(h, :), 'LineWidth', 1.5, ...
			'DisplayName', 'Dynamical Dependence', 'LineStyle', '-', 'Color', colors(color_idx, :));
		hold on;
		color_idx = color_idx + 1;
	end
	ylabel('Dynamical Dependence', 'FontSize', 15, 'interpreter', 'latex');
	set(gca, 'YColor', 'k');

	% Add margin to right y-axis
	ylim_right = ylim;
	margin_right = 0.1 * (ylim_right(2) - ylim_right(1));
	ylim([ylim_right(1) - margin_right, ylim_right(2) + margin_right]);

	% Common plot settings
	xlim([min(all_phase_lag), max(all_phase_lag)]);
	xlabel('Phase Lag', 'FontSize', 15, 'interpreter', 'latex');

	% Create title with explicit handle and properties
	title({[config.title ' vs Phase Lag'], ...
		['RMI = ' num2str(all_rmi(end)) ', Coupling Ratio Mean = ' num2str(actual_coup_ratio, '%.1f')]}, ...
		'interpreter', 'latex');

	legend('Location', [0.2, 0.15, 0.2, 0.1], 'interpreter', 'latex', 'FontSize', 12);
	grid on;
	set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

	% Force title font size update
	ax = gca;
	ax.Title.FontSize = 20;

	% Save plot
	filename_lineplot = [network '_lineplot_coup_ratio_' ...
		strrep(num2str(actual_coup_ratio, '%.1f'), '.', 'p') '_DD_' config.filename_suffix];
	location = string(strcat(pathout_plots_measures, filename_lineplot, '.png'));
	exportgraphics(gcf, location);
end