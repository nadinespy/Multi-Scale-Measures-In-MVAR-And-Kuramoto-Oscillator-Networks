
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_chi_met_km_rand_coup_ratio_phase_lag_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In order to run this script, 
% - mec_km_rand_coup_ratio_phase_lag_rmi_setup.m
% - params_chi_met_rand_coup_ratio_phase_lag_rmi.m 
% (or any other script with the prefix 'params')
% must be run (in that order).
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define model name
network = [num2str(n_nodes) 'km' '_lag' num2str(ar_model_order)];

% filename for this particular model & measure
filename = strcat(num2str(n_nodes), 'km_lag', num2str(ar_model_order), ...
	'_rand_coup_ratio_phase_lag_rmi_', 'chi_met');
load([char(pathout_data_measures), char(filename) '.mat'], ...
	'results');

% Nasty work-around: if RMI dimension entails only one value, 
% 'results' will not be n_phase_lag x n_coup_ratio_mean x n_rmi, 
% but just n_phase_lag x n_coup_ratio_mean (singleton dimensions 
% past the 2nd dimension are automatically truncated), hence the 
% need to manually add a first singleton dimension to make the 
% below loop for plotting work.
if ndims(results) == 2
	new_results(1,:,:) = results;
else 
	new_results = results;
end

measure_in_struct = {'global_metastability', 'global_chimera_index'};
measure_in_title = {'Metastability', 'Chimera Index'};

%% HEATMAP FOR AVERAGE MEASURE (ONE FOR EACH COUPLING RATIO MEAN) 

% custom function to create new axes labels where the xtick label 
% is displayed only every number of times (number is given by 
% phase_lag_ticks_steps/coup_ratio_mean_ticks_steps) 
new_phase_lag_axis = create_axis_ticks(all_phase_lag, ...
	phase_lag_ticks_steps);
new_coup_ratio_mean_axis = create_axis_ticks(all_coup_ratio_mean, ...
	coup_ratio_mean_ticks_steps);

for g = 1:length(measure_in_struct)

	% get means of DDCommSync
	measure_means = [];
	for h = 1:n_rmi_to_plot
		for i = 1:n_coup_ratio_mean
			for j = 1:n_phase_lag
				measure_means(end+1) = ...
					mean(mean(new_results(h,i,j). ...
					(measure_in_struct{g})));
			end
		end
	end
	% find global min/max across all slices
	clim_min = min(measure_means);
	clim_max = max(measure_means);

	figure;
	for h = 1:n_rmi_to_plot
    	subplot(1, n_rmi_to_plot, h)
    	measure_slice = zeros(n_samples_noise_corr, n_samples_noise_corr);
    	
    	for i = 1:n_coup_ratio_mean
        	for j = 1:n_phase_lag
            	measure_slice(i,j) = mean(mean(results(i,j).(measure_in_struct{g})));
        	end
	end
    	
	imagesc(measure_slice);
    	colormap(gca, parula);
	set(gca, 'CLim', [clim_min clim_max]);

	colorbar;

    	% rest of the plotting settings (no more CLim changes)
    	yticks(1:length(new_coup_ratio_mean_axis));
    	xticks(1:length(new_phase_lag_axis));
    	set(gca, 'TickLength', [0 0]);
    	xticklabels(new_phase_lag_axis);
    	yticklabels(new_coup_ratio_mean_axis);
    	xtickangle(45);
    	ytickangle(45);
    	ax = gca;
    	ax.YDir = 'normal';

	% adjust font sizes
	set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');  % for tick labels
	title(['RMI = ' num2str(all_rmi(h))], 'FontSize', 17, 'interpreter','latex');  
	ylabel('coupling ratio mean', 'FontSize', 15, 'interpreter','latex');  
	xlabel('phase lag', 'FontSize', 15, 'interpreter','latex');  

	cb = colorbar;
	set(cb, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

	end

	% sgtitle([measure_in_title], 'FontSize', 14);
	sgtitle(measure_in_title{g}, 'FontSize', 19, 'interpreter','latex');

	filename_heatmap = {[network ...
		'_heatmap_coup_ratio_phase_lag_rmi_', ...
		measure_in_struct{g}]};
	
	location = string(strcat(pathout_plots_measures, ...
		filename_heatmap, '.png'));
	exportgraphics(gcf, location);
end

clear results new_results;

%% LINE PLOTS - METASTABILITY/CHIMERA INDEX & DD FOR SPECIFIC COUPLING RATIO MEANS

% Define target coupling ratio means
target_coup_ratios = [0.2, 0.9];

% Find the indices corresponding to target coupling ratio means
coup_ratio_indices = zeros(size(target_coup_ratios));
actual_coup_ratios = zeros(size(target_coup_ratios));
for i = 1:length(target_coup_ratios)
    [~, coup_ratio_indices(i)] = min(abs(all_coup_ratio_mean - target_coup_ratios(i)));
    actual_coup_ratios(i) = all_coup_ratio_mean(coup_ratio_indices(i));
end

% Loop through each coupling ratio
for coup_idx = 1:length(target_coup_ratios)
    
    % ------------------------------------------------
    % LOAD AND PLOT CHI_MET DATA
    % ------------------------------------------------

    fprintf('Loading chi_met data for coupling ratio %.1f...\n', target_coup_ratios(coup_idx));
    
    % Load chi_met data
    filename_chi_met = strcat(num2str(n_nodes), 'km_lag', num2str(ar_model_order), ...
        '_rand_coup_ratio_phase_lag_rmi_', 'chi_met');
    load([char(pathout_data_measures), char(filename_chi_met) '.mat'], 'results');
    
    % Handle singleton dimension issue
    if ndims(results) == 2
        chi_met_results(1,:,:) = results;
    else 
        chi_met_results = results;
    end
    
    % Extract data for both measures
    chi_met_measures = {'global_metastability', 'global_chimera_index'};
    chi_met_titles = {'Metastability', 'Chimera Index'};
    chi_met_data = cell(1, 2); % Store data for both measures
    
    for g = 1:length(chi_met_measures)
        chi_met_data{g} = zeros(n_rmi_to_plot, n_phase_lag);
        for h = 1:n_rmi_to_plot
            for j = 1:n_phase_lag
                chi_met_data{g}(h, j) = mean(mean(chi_met_results(h, coup_ratio_indices(coup_idx), j).(chi_met_measures{g})));
            end
        end
    end
    
    clear results chi_met_results; % Free memory
    
    % ------------------------------------------------
    % LOAD DD DATA
    % ------------------------------------------------

    fprintf('Loading DD data for coupling ratio %.1f...\n', target_coup_ratios(coup_idx));
    
    % Load DD data
    filename_dd = strcat(num2str(n_nodes), 'km_lag', num2str(ar_model_order), ...
        '_rand_coup_ratio_phase_lag_rmi_', 'dd_ce_co_info');
    load([char(pathout_data_measures), char(filename_dd) '.mat'], 'results');
    
    % Handle singleton dimension issue
    if ndims(results) == 2
        dd_results(1,:,:) = results;
    else 
        dd_results = results;
    end
    
    % Extract DD data
    dd_data = zeros(n_rmi_to_plot, n_phase_lag);
    for h = 1:n_rmi_to_plot
        for j = 1:n_phase_lag
            dd_data(h, j) = mean(mean(dd_results(h, coup_ratio_indices(coup_idx), j).DDCommSync));
        end
    end
    
    clear results dd_results; % Free memory
    
    % ------------------------------------------------
    % CREATE PLOTS
    % ------------------------------------------------


    % Plot 1: DD and Chimera Index (dual y-axes)
    figure;
    colors = lines(n_rmi_to_plot * 2); % Get enough colors for all lines
    color_idx = 1;

    % Create first y-axis for Chimera Index
    yyaxis left;
    for h = 1:n_rmi_to_plot
	    plot(all_phase_lag, chi_met_data{2}(h, :), 'LineWidth', 1.5, ...
		    'DisplayName', 'Chimera Index', ...
		    'LineStyle', '-', 'Color', colors(color_idx, :));
	    hold on;
	    color_idx = color_idx + 1;
    end
    ylabel('Chimera Index', 'FontSize', 15, 'interpreter', 'latex');
    set(gca, 'YColor', 'k'); % Keep left axis black

    % Add margin to left y-axis
    ylim_left = ylim;
    margin_left = 0.1 * (ylim_left(2) - ylim_left(1));
    ylim([ylim_left(1) - margin_left, ylim_left(2) + margin_left]);

    % Create second y-axis for DD
    yyaxis right;
    for h = 1:n_rmi_to_plot
	    plot(all_phase_lag, dd_data(h, :), 'LineWidth', 1.5, ...
		    'DisplayName', 'Dynamical Dependence', ...
		    'LineStyle', '-', 'Color', colors(color_idx, :));
	    hold on;
	    color_idx = color_idx + 1;
    end
    ylabel('Dynamical Dependence', 'FontSize', 15, 'interpreter', 'latex');
    set(gca, 'YColor', 'k'); % Keep right axis black

    % Add margin to right y-axis
    ylim_right = ylim;
    margin_right = 0.1 * (ylim_right(2) - ylim_right(1));
    ylim([ylim_right(1) - margin_right, ylim_right(2) + margin_right]);

    xlim([min(all_phase_lag), max(all_phase_lag)]);
    xlabel('Phase Lag', 'FontSize', 15, 'interpreter', 'latex');
    title(['DD and Chimera Index vs Phase Lag, RMI = ' num2str(all_rmi(h)) ', Coupling Ratio Mean = ' num2str(actual_coup_ratios(coup_idx), '%.1f')], ...
	    'FontSize', 17, 'interpreter', 'latex');

    legend('Location', [0.2, 0.15, 0.2, 0.1], 'interpreter', 'latex', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    % Save plot
    filename_lineplot = [network '_lineplot_coup_ratio_' strrep(num2str(actual_coup_ratios(coup_idx), '%.1f'), '.', 'p') '_DD_ChimeraIndex'];
    location = string(strcat(pathout_plots_measures, filename_lineplot, '.png'));
    exportgraphics(gcf, location);

    % Plot 2: DD and Metastability (dual y-axes)
    figure;
    colors = lines(n_rmi_to_plot * 2); % Get enough colors for all lines
    color_idx = 1;

    % Create first y-axis for Metastability
    yyaxis left;
    for h = 1:n_rmi_to_plot
	    plot(all_phase_lag, chi_met_data{1}(h, :), 'LineWidth', 1.5, ...
		    'DisplayName', 'Metastability', ...
		    'LineStyle', '-', 'Color', colors(color_idx, :));
	    hold on;
	    color_idx = color_idx + 1;
    end
    ylabel('Metastability', 'FontSize', 15, 'interpreter', 'latex');
    set(gca, 'YColor', 'k'); % Keep left axis black

    % Add margin to left y-axis
    ylim_left = ylim;
    margin_left = 0.1 * (ylim_left(2) - ylim_left(1));
    ylim([ylim_left(1) - margin_left, ylim_left(2) + margin_left]);

    % Create second y-axis for DD
    yyaxis right;
    for h = 1:n_rmi_to_plot
	    plot(all_phase_lag, dd_data(h, :), 'LineWidth', 1.5, ...
		    'DisplayName', 'Dynamical Dependence', ...
		    'LineStyle', '-', 'Color', colors(color_idx, :));
	    hold on;
	    color_idx = color_idx + 1;
    end
    ylabel('Dynamical Dependence', 'FontSize', 15, 'interpreter', 'latex');
    set(gca, 'YColor', 'k'); % Keep right axis black

    % Add margin to right y-axis
    ylim_right = ylim;
    margin_right = 0.1 * (ylim_right(2) - ylim_right(1));
    ylim([ylim_right(1) - margin_right, ylim_right(2) + margin_right]);

    xlim([min(all_phase_lag), max(all_phase_lag)]);
    xlabel('Phase Lag', 'FontSize', 15, 'interpreter', 'latex');
    title(['DD and Metastability vs Phase Lag, RMI = ' num2str(all_rmi(h)) ', Coupling Ratio Mean = ' num2str(actual_coup_ratios(coup_idx), '%.1f')], ...
	    'FontSize', 17, 'interpreter', 'latex');

    legend('Location', [0.2, 0.15, 0.2, 0.1], 'interpreter', 'latex', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    % Save plot
    filename_lineplot = [network '_lineplot_coup_ratio_' strrep(num2str(actual_coup_ratios(coup_idx), '%.1f'), '.', 'p') '_DD_Metastability'];
    location = string(strcat(pathout_plots_measures, filename_lineplot, '.png'));
    exportgraphics(gcf, location);

    % Clear data for this coupling ratio before moving to next
    clear chi_met_data dd_data;

    fprintf('Completed plots for coupling ratio %.1f\n', actual_coup_ratios(coup_idx));
end

fprintf('All line plots completed!\n');
fprintf('Actual coupling ratios used: %.3f and %.3f\n', actual_coup_ratios(1), actual_coup_ratios(2));




%% LINE PLOTS - METASTABILITY/CHIMERA INDEX & DD FOR SPECIFIC COUPLING RATIO MEANS

% Define target coupling ratio means
target_coup_ratios = [0.2, 0.9];

% Find the indices corresponding to target coupling ratio means
coup_ratio_indices = zeros(size(target_coup_ratios));
actual_coup_ratios = zeros(size(target_coup_ratios));
for i = 1:length(target_coup_ratios)
    [~, coup_ratio_indices(i)] = min(abs(all_coup_ratio_mean - target_coup_ratios(i)));
    actual_coup_ratios(i) = all_coup_ratio_mean(coup_ratio_indices(i));
end

% Define plot configurations
plot_configs = {
    struct('data_idx', 2, 'name', 'Chimera Index', 'title', 'DD and Chimera Index', 'filename_suffix', 'ChimeraIndex'), ...
    struct('data_idx', 1, 'name', 'Metastability', 'title', 'DD and Metastability', 'filename_suffix', 'Metastability')
};

% Loop through each coupling ratio
for coup_idx = 1:length(target_coup_ratios)
    
    % ------------------------------------------------
    % LOAD AND EXTRACT CHI_MET DATA
    % ------------------------------------------------
    
    fprintf('Loading chi_met data for coupling ratio %.1f...\n', target_coup_ratios(coup_idx));
    
    % Load chi_met data
    filename_chi_met = strcat(num2str(n_nodes), 'km_lag', num2str(ar_model_order), ...
        '_rand_coup_ratio_phase_lag_rmi_', 'chi_met');
    load([char(pathout_data_measures), char(filename_chi_met) '.mat'], 'results');
    
    % Handle singleton dimension issue
    if ndims(results) == 2
        chi_met_results(1,:,:) = results;
    else 
        chi_met_results = results;
    end
    
    % Extract data for both measures
    chi_met_measures = {'global_metastability', 'global_chimera_index'};
    chi_met_data = cell(1, 2);
    
    for g = 1:length(chi_met_measures)
        chi_met_data{g} = zeros(n_rmi_to_plot, n_phase_lag);
        for h = 1:n_rmi_to_plot
            for j = 1:n_phase_lag
                chi_met_data{g}(h, j) = mean(mean(chi_met_results(h, coup_ratio_indices(coup_idx), j).(chi_met_measures{g})));
            end
        end
    end
    
    clear results chi_met_results;
    
    % ------------------------------------------------
    % LOAD AND EXTRACT DD DATA
    % ------------------------------------------------

    fprintf('Loading DD data for coupling ratio %.1f...\n', target_coup_ratios(coup_idx));
    
    % Load DD data
    filename_dd = strcat(num2str(n_nodes), 'km_lag', num2str(ar_model_order), ...
        '_rand_coup_ratio_phase_lag_rmi_', 'dd_ce_co_info');
    load([char(pathout_data_measures), char(filename_dd) '.mat'], 'results');
    
    % Handle singleton dimension issue
    if ndims(results) == 2
        dd_results(1,:,:) = results;
    else 
        dd_results = results;
    end
    
    % Extract DD data
    dd_data = zeros(n_rmi_to_plot, n_phase_lag);
    for h = 1:n_rmi_to_plot
        for j = 1:n_phase_lag
            dd_data(h, j) = mean(mean(dd_results(h, coup_ratio_indices(coup_idx), j).DDCommSync));
        end
    end
    
    clear results dd_results;
    
    % ------------------------------------------------
    % CREATE PLOTS USING FUNCTION
    % ------------------------------------------------
    
    % Create both plots using the function
    for plot_idx = 1:length(plot_configs)
        plot_dual_lines(all_phase_lag, chi_met_data, dd_data, plot_configs{plot_idx}, ...
            actual_coup_ratios(coup_idx), all_rmi, n_rmi_to_plot, network, pathout_plots_measures);
    end
    
    % Clear data for this coupling ratio before moving to next
    clear chi_met_data dd_data;
    
    fprintf('Completed plots for coupling ratio %.1f\n', actual_coup_ratios(coup_idx));
end

fprintf('All line plots completed!\n');
fprintf('Actual coupling ratios used: %.3f and %.3f\n', actual_coup_ratios(1), actual_coup_ratios(2));
