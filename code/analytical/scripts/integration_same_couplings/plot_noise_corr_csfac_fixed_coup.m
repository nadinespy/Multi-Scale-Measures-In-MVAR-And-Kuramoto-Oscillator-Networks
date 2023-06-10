%% PLOTTING: INTEGRATION IN TWO-NODE MVAR NETWORKS - NOISE CORRELATION & C-FACTOR, FIXED COUPLING

%% line plots
% -------------------------------------------------------------------------
% {
% plot INTEGRATION (y-axis) and C-FACTOR (x-axis), one line per NOISE CORRELATION, FIXED COUPLING
figure(1); clf
plot(noise_corrs',I');

lgd = legend(num2str(noise_corrs'));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs c-scaling (n = %d, model order = %d, coupling = %d),\n single line for each noise correlation',n,p,two_node_couplings(coup_index)));
xlabel('connectivity scale factor');
ylabel('integration (nats)');

grid on

file_names = {['line_noise_corr_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);

% plot INTEGRATION (y-axis) and NOISE CORRELATION (x-axis), one line per C-FACTORS, FIXED COUPLING
figure(1); clf
plot(noise_corrs',I);

lgd = legend(num2str(csfac));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs noise correlation (n = %d, model order = %d, coupling = %d),\n single line for each c-factor',n,p,two_node_couplings(coup_index)));
xlabel('noise correlation');
ylabel('integration (nats)');

grid on

file_names = {['line_csfac_noise_corr_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);
%}

%% heatmaps

% {
% -------------------------------------------------------------------------
y_axis_noise_corr = {'0.0', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.32', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.56', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.9'};

x_axis_csfac = {'0.0', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.33', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.66', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.99'};

% get column names for table
c_factor_str = {};
for t = 1:length(csfac)
	c_factor_str{t} = num2str(csfac(t));
end

% get row names for table
noise_corr_str = {};
for t = 1:length(noise_corrs)
	noise_corr_str{t} = num2str(noise_corrs(t));
end

y_label_noise_corr	= 'noise correlation';
x_label_csfac		= 'c-factor';

% plot INTEGRATION for all values NOISE CORRELATION (y-axis) and C-FACTORS (x-axis), FIXED COUPLING
file_names_heatmap = {['heatmap_int_noise_corr_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('integration: noise correlation vs c-factor (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};

% get row names for table
	noise_corr_str = {};
	for t = 1:length(noise_corrs)
		noise_corr_str{t} = num2str(noise_corrs(t));
	end
	
integration = array2table(I, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
integration = {integration};

plot_heatmap(integration, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);

% plot DOUBLE-RED MMI for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_rtr_mmi_noise_corr_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-red MMI: noise correlation vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidRed_MMI_plot = array2table(phiidRed_MMI, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
phiidRed_MMI_plot = {phiidRed_MMI_plot};

plot_heatmap(phiidRed_MMI_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);

% plot DOUBLE-RED CCS for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_rtr_ccs_noise_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-red CCS: noise correlation vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidRed_CCS_plot = array2table(phiidRed_CCS, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
phiidRed_CCS_plot = {phiidRed_CCS_plot};

plot_heatmap(phiidRed_CCS_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);

% plot CE MMI for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_ce_mmi_noise_corr_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('CE MMI: noise correlation vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidCE_MMI_plot = array2table(phiidCE_MMI, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
phiidCE_MMI_plot = {phiidCE_MMI_plot};

plot_heatmap(phiidCE_MMI_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);

% plot CE CCS for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_ce_ccs_noise_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('CE CCS: noise correlation vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidCE_CCS_plot = array2table(phiidCE_CCS, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
phiidCE_CCS_plot = {phiidCE_CCS_plot};

plot_heatmap(phiidCE_CCS_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);

% plot DOUBLE-SYN MMI for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_sts_mmi_noise_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-syn MMI: noise correlation vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidSyn_MMI_plot = array2table(phiidSyn_MMI, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
phiidSyn_MMI_plot = {phiidSyn_MMI_plot};

plot_heatmap(phiidSyn_MMI_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);

% plot DOUBLE-SYN CCS for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_sts_ccs_noise_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-syn CCS: noise correlation vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidSyn_CCS_plot = array2table(phiidSyn_CCS, 'RowNames', noise_corr_str,'VariableNames', c_factor_str);
phiidSyn_CCS_plot = {phiidSyn_CCS_plot};

plot_heatmap(phiidSyn_CCS_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_noise_corr, x_label_csfac, y_label_noise_corr, pathout_plots);
%}