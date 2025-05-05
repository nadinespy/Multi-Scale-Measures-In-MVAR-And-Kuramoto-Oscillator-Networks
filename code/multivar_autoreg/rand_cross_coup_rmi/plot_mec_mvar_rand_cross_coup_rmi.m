
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_mec_mvar_rand_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is run via plot_all_mec_mvar_rand_cross_coup_rmi.m.
% The following scripts must be run beforehand:
% - params_phiid_mmi_rand_cross_coup_rmi.m 
% (or any other script with the prefix 'params'),
% - mec_mvar_rand_cross_coup_rmi_local_dirs.m
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define model name
network = [num2str(n_nodes) 'mvar' '_lag' num2str(ar_model_order)];

% filename for this particular model & measure
filename = strcat(num2str(n_nodes), 'mvar_lag', num2str(ar_model_order), ...
	'_rand_cross_coup_rmi_', measure_in_filename);
load([char(pathout_data_measures), char(filename) '.mat'], ...
	'results');

% define measure name as in table
if strcmp(measure_in_filename, 'phiid_mmi_ce') 
	measure_in_struct = {'PhiID_CE_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Emergence Capacity using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_dc') 
	measure_in_struct = {'PhiID_DC_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Downward Causation using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_cd') 
	measure_in_struct = {'PhiID_CD_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Causal Decoupling using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_uc') 
	measure_in_struct = {'PhiID_UC_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Upward Causation using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_double_red') 
	measure_in_struct = {'PhiID_DoubleRed_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Double-Redundancy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_double_syn') 
	measure_in_struct = {'PhiID_DoubleSyn_MMI'};
	measure_in_title = {'\$\Phi\textrm{ID}$-based Double-Redundancy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_syn') 
	measure_in_struct = {'PhiID_Syn_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Synergy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_transfer') 
	measure_in_struct = {'PhiID_Transfer_MMI'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Transfer Entropy using MMI'};
elseif strcmp(measure_in_filename, 'phiid_ccs_ce') 
	measure_in_struct = {'PhiID_CE_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Emergence Capacity using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_dc') 
	measure_in_struct = {'PhiID_DC_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Downward Causation using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_cd') 
	measure_in_struct = {'PhiID_CD_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Causal Decoupling using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_uc') 
	measure_in_struct = {'PhiID_UC_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Upward Causation using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_double_red') 
	measure_in_struct = {'PhiID_DoubleRed_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Double-Redundancy using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_double_syn')
	measure_in_struct = {'PhiID_DoubleSyn_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Synergy using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_syn') 
	measure_in_struct = {'PhiID_Syn_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Synergy using CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_transfer') 
	measure_in_struct = {'PhiID_Transfer_CCS'};
	measure_in_title = {'$\Phi\textrm{ID}$-based Transfer Entropy using CCS'};
elseif strcmp(measure_in_filename, 'dd_ce_co_info')
	measure_in_struct = {'DD_PCA', 'ShannonCE_PCA', ...
		'DDGrassMin', 'DDGrassMean', 'ShannonCEGrassMax', ...
		'ShannonCEGrassMean', 'CoInfoPCA', 'CoInfoGrassMin', ...
		'CoInfoGrassMax', 'CoInfoGrassMean'};
	measure_in_title = {{'Dynamical Dependence', 'for 1st Principal Component'}, ...
		{'Whole-Part-Emergence for', '1st Principal Component'}, ...
		{'Dynamical Dependence for', 'DD-optimised linear coarse-graining'}, ...
		{'Dynamical Dependence for', 'DD-optimised linear coarse-graining'}, ...
		{'Whole-Part-Emergence for', ...
		'WPE-optimised linear coarse-graining'}, ...
		{'Whole-Part-Emergence for', ...
		'WPE-optimised linear coarse-graining'}, ...
		{'Co-Information for 1st Principal Component'}, ...
		{'Co-Information for CI-optimised', 'linear coarse-graining'}, ...
		{'Co-Information for CI-optimised', 'linear coarse-graining'}, ...
		{'Co-Information for CI-optimised', 'linear coarse-graining'}};
elseif strcmp(measure_in_filename, 'integrated_info')
	measure_in_struct = {'IntegratedInfo', ...
		'CausalDensity', ...
		'IntegratedSynergy', ...
		'DecoderIntegration', ...
		'IntegratedInteraction'};
	measure_in_title = {'Integrated Information', ...
		'Causal Density', 'Integrated Synergy', ...
		'Decoder-Based Integrated Information', ...
		'Integrated Interaction'};
elseif strcmp(measure_in_filename, 'control')
	measure_in_struct = {'TimeDelayedMI'}; %{'TimeDelayedMI', 'MultiInfo', 'AverageCorr'};
	measure_in_title = {'Time-Delayed Mutual Information'}; 
		%{'Time-Delayed Mutual Information', ...
		%'Multi-Information', ...
		%'Average Correlation'};
else 
	error('Unknown measure ''%s''', measure_in_filename);
end

%% HEATMAP FOR AVERAGE MEASURE (ONE FOR EACH COUPLING RATIO MEAN) 

% custom function to create new axes labels where the xtick label 
% is displayed only every number of times (number is given by 
% phase_lag_ticks_steps/rmi_ticks_steps) 
new_cross_coup_axis = create_axis_ticks(all_cross_coup, ...
	cross_coup_ticks_steps);
new_rmi_axis = create_axis_ticks(all_rmi, rmi_ticks_steps);

% create a figure with subplots for each coupling ratio mean value
for g = 1:length(measure_in_struct)

	% get means of measure
	for i = 1:n_cross_coup
		for j = 1:n_rmi
			measure(i,j) = mean(mean(results(i,j).(measure_in_struct{g})));
			all_global_coup(i,j) = mean(mean(results(i,j).global_coup));
		end
	end

	% create evenly spaced vector within range of global coupling
	% magnitudes
	min_global_coup = min(min(all_global_coup));
	max_global_coup = max(max(all_global_coup));
	all_global_coup_for_plot = round(linspace(min_global_coup, ...
		max_global_coup, length(all_global_coup)),2);

	% plot
	figure;

	imagesc(measure);
    	colormap(gca, parula);
	clim = [min(measure(:)), max(measure(:))];
	set(gca, 'CLim', clim);

	colorbar;
	
	xticks(1:length(new_rmi_axis));
	yticks(1:length(new_cross_coup_axis));
	set(gca, 'TickLength', [0 0]);
	yticklabels(new_cross_coup_axis);
	xticklabels(new_rmi_axis);
	xtickangle(45);
	ytickangle(45);

	ax = gca;
	ax.YDir = 'normal';

	% adjust font sizes
	set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');  % for tick labels
	cb = colorbar;
	set(cb, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
		
	title(measure_in_title{g}, 'FontSize', 17, 'interpreter','latex');
	xlabel('RMI', 'FontSize', 15, 'interpreter','latex');
	ylabel('global coupling', 'FontSize', 15, 'interpreter','latex');

	filename_heatmap = {[network '_heatmap_cross_coup_rmi_', ...
		measure_in_struct{g}]};

	location = string(strcat(pathout_plots_measures, filename_heatmap, '.png'));
	exportgraphics(gcf, location);
	
	close all;

	% {
	if g == length(measure_in_struct)
		clear results;
	end
	%}
end
