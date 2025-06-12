
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_mec_km_rand_coup_ratio_phase_lag_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In order to run this script, 
% - params_phiid_mmi_rand_coup_ratio_phase_lag_rmi.m 
% (or any other script with the prefix 'params'),
% - mec_km_rand_coup_ratio_phase_lag_rmi_local_dirs.m
% must be run (in that order).
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define model name
network = [num2str(n_nodes) 'km' '_lag' num2str(ar_model_order)];

% filename for this particular model & measure
filename = strcat(num2str(n_nodes), 'km_lag', num2str(ar_model_order), ...
	'_rand_coup_ratio_phase_lag_rmi_', measure_in_filename);
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

% define measure name as in table
if strcmp(measure_in_filename, 'multi_info') 
	measure_in_struct = {'MultiInfo'};
	measure_in_title = {'Multi-Information'};
elseif strcmp(measure_in_filename, 'phiid_mmi_ce') 
	measure_in_struct = {'PhiID_CE_MMI'};
	measure_in_title = {'\PhiID-based Whole-Parts-Emergence using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_dc') 
	measure_in_struct = {'PhiID_DC_MMI'};
	measure_in_title = {'\PhiID-based Downward Causation using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_cd') 
	measure_in_struct = {'PhiID_CD_MMI'};
	measure_in_title = {'\PhiID-based Causal Decoupling using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_uc') 
	measure_in_struct = {'PhiID_UC_MMI'};
	measure_in_title = {'\PhiID-based Upward Causation using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_double_red') 
	measure_in_struct = {'PhiID_DoubleRed_MMI'};
	measure_in_title = {'\PhiID-based Double-Redundancy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_double_syn') 
	measure_in_struct = {'PhiID_DoubleSyn_MMI'};
	measure_in_title = {'\PhiID-based Double-Redundancy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_syn') 
	measure_in_struct = {'PhiID_Syn_MMI'};
	measure_in_title = {'\PhiID-based Synergy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_transfer') 
	measure_in_struct = {'PhiID_Transfer_MMI'};
	measure_in_title = {'\PhiID-based Transfer Entropy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_ccs_ce') 
	measure_in_struct = {'PhiID_CE_CCS'};
	measure_in_title = {'\PhiID-based Whole-Parts-Emergence using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_dc') 
	measure_in_struct = {'PhiID_DC_CCS'};
	measure_in_title = {'\PhiID-based Downward Causation using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_cd') 
	measure_in_struct = {'PhiID_CD_CCS'};
	measure_in_title = {'\PhiID-based Causal Decoupling using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_uc') 
	measure_in_struct = {'PhiID_UC_CCS'};
	measure_in_title = {'\PhiID-based Upward Causation using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_double_red') 
	measure_in_struct = {'PhiID_DoubleRed_CCS'};
	measure_in_title = {'\PhiID-based Double-Redundancy using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_double_syn')
	measure_in_struct = {'PhiID_DoubleSyn_CCS'};
	measure_in_title = {'\PhiID-based Synergy using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_syn') 
	measure_in_struct = {'PhiID_Syn_CCS'};
	measure_in_title = {'\PhiID-based Synergy using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_transfer') 
	measure_in_struct = {'PhiID_Transfer_CCS'};
	measure_in_title = {'\PhiID-based Transfer Entropy using CCS'};
elseif strcmp(measure_in_filename, 'dd_ce_co_info')
	measure_in_struct = {'DDCommSync'};
	measure_in_title = {{'Dynamical Dependence for Community Synchronies'}};
elseif strcmp(measure_in_filename, 'integrated_info')
	measure_in_struct = {'IntegratedInfo', ...
		'CausalDensity', ...
		'IntegratedSynergy', ...
		'DecoderIntegration', ...
		'IntegratedInteraction'};
	measure_in_title = {'Integrated Information', 'Causal Density', 'Integrated Synergy', ...
		'Decoder-Based Integrated Information', 'Integrated Interaction'};
elseif strcmp(measure_in_filename, 'control')
	measure_in_struct = {'TimeDelayedMI', 'AverageCorr'};
	measure_in_title = {'Time-Delayed Mutual Information', 'Average Correlation'};
else 
	error('Unknown measure ''%s''', measure_in_filename);
end

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

%% LINE PLOT FOR COUPLING RATIO MEAN = 0.2
% Find the index corresponding to coupling ratio mean = 0.2
target_coup_ratio_mean = 0.2;
[~, coup_ratio_idx] = min(abs(all_coup_ratio_mean - target_coup_ratio_mean));
actual_coup_ratio_mean = all_coup_ratio_mean(coup_ratio_idx);

% Create line plots for each measure
for g = 1:length(measure_in_struct)
    figure;
    
    % Plot for each RMI value
    for h = 1:n_rmi_to_plot
        % Extract data for the specific coupling ratio mean across all phase lags
        line_data = zeros(1, n_phase_lag);
        for j = 1:n_phase_lag
            line_data(j) = mean(mean(new_results(h, coup_ratio_idx, j).(measure_in_struct{g})));
        end
        
        % Plot the line
        plot(all_phase_lag, line_data, 'LineWidth', 1.5, 'DisplayName', ['RMI = ' num2str(all_rmi(h))]);
        hold on;
    end
    
    % Set x-axis limits to match the data range exactly
    xlim([min(all_phase_lag), max(all_phase_lag)]);
    
    % Formatting
    xlabel('Phase Lag', 'FontSize', 15, 'interpreter', 'latex');
    ylabel(measure_in_title{g}, 'FontSize', 15, 'interpreter', 'latex');
    title([measure_in_title{g} ' vs Phase Lag (Coupling Ratio Mean = ' num2str(actual_coup_ratio_mean, '%.1f') ')'], ...
        'FontSize', 17, 'interpreter', 'latex');
    
    % Add legend if plotting multiple RMI values
    if n_rmi_to_plot > 1
        legend('Location', 'best', 'interpreter', 'latex', 'FontSize', 12);
    end
    
    grid on;
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
    
    % Save the plot
    filename_lineplot = [network '_lineplot_coup_ratio_0p2_phase_lag_' measure_in_struct{g}];
    location = string(strcat(pathout_plots_measures, filename_lineplot, '.png'));
    exportgraphics(gcf, location);
end

% Display the actual coupling ratio mean value used
fprintf('Line plot created for coupling ratio mean = %.3f (closest to requested 0.2)\n', actual_coup_ratio_mean);

clear results;
