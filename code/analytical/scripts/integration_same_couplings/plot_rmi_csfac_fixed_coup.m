%% PLOTTING: INTEGRATION IN TWO-NODE MVAR NETWORKS - RMI, FIXED COUPLING

%% line plots

%{
% -------------------------------------------------------------------------
% plot INTEGRATION and C-FACTOR, and FIXED RMI and FIXED COUPLING
% with corresponding total and generalized variance
figure(1); clf
plot(csfac,[I(rmi_index,:) Iv(rmi_index,:) Ig(rmi_index,:)]);
yline(all_rmi(rmi_index),'g','residuals multi-information');
title(sprintf('Integration vs connectivity scaling (n = %d, model order = %d, rmi = %g, coupling = %d)\n',n,p,all_rmi(rmi_index), two_node_couplings(coup_index)));
legend({'integration','total variance','generalised variance'});
xlabel('connectivity scale factor');
ylabel('integration (nats)');

grid on

file_names = {['line_2node_csfac_rmi_' num2str(all_rmi(rmi_index)) '_coup_' num2str(two_node_couplings(coup_index))]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);
%}

% {
% -------------------------------------------------------------------------
% plot INTEGRATION (y-axis) and C-FACTOR (x-axis), one line per RMI, FIXED COUPLING
figure(1); clf
plot(all_rmi',I'); % plots single line for each row 

lgd = legend(num2str(all_rmi'));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs c-scaling (n = %d, model order = %d, coupling = %d),\n single line for each rmi', ...
	n,time_lag,two_node_couplings(coup_index)));
xlabel('connectivity scale factor');
ylabel('integration (nats)');

grid on

file_names = {['line_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);

% plot INTEGRATION (y-axis) and RMI (x-axis), one line per C-FACTOR, FIXED COUPLING
figure(1); clf
plot(all_rmi',I);

lgd = legend(num2str(csfac));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs rmi (n = %d, model order = %d, coupling = %d),\n single line for each c-factor',n,p,two_node_couplings(coup_index)));
xlabel('rmi');
ylabel('integration (nats)');

grid on

file_names = {['line_csfac_rmi_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);
%}

%% heatmaps

% {
% -------------------------------------------------------------------------

y_axis_rmi = {'0.0', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.33', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.66', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '1'};

x_axis_csfac = {'0.0', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.33', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.66', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.99'};

% get row names for table
all_rmi_str = {};
for e = 1:length(all_rmi)
	all_rmi_str{e} = num2str(all_rmi(e));
end

% get column names for table
c_factor_str = {};
for t = 1:length(csfac)
	c_factor_str{t} = num2str(csfac(t));
end

y_label_rmi			= 'rmi';
x_label_csfac		= 'c-factor';

% plot INTEGRATION for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_int_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('integration: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
I_plot = array2table(I, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
I_plot = {I_plot};

plot_heatmap(I_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);

% plot DOUBLE-RED MMI for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_rtr_mmi_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-red MMI: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidRed_MMI_plot = array2table(phiidRed_MMI, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
phiidRed_MMI_plot = {phiidRed_MMI_plot};

plot_heatmap(phiidRed_MMI_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);

% plot DOUBLE-RED CCS for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_rtr_ccs_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-red CCS: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidRed_CCS_plot = array2table(phiidRed_CCS, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
phiidRed_CCS_plot = {phiidRed_CCS_plot};

plot_heatmap(phiidRed_CCS_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);

% plot CE MMI for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_ce_mmi_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('CE MMI: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidCE_MMI_plot = array2table(phiidCE_MMI, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
phiidCE_MMI_plot = {phiidCE_MMI_plot};

plot_heatmap(phiidCE_MMI_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);

% plot CE CCS for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_ce_ccs_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('CE CCS: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidCE_CCS_plot = array2table(phiidCE_CCS, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
phiidCE_CCS_plot = {phiidCE_CCS_plot};

plot_heatmap(phiidCE_CCS_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);

% plot DOUBLE-SYN MMI for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_sts_mmi_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-syn MMI: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidSyn_MMI_plot = array2table(phiidSyn_MMI, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
phiidSyn_MMI_plot = {phiidSyn_MMI_plot};

plot_heatmap(phiidSyn_MMI_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);

% plot DOUBLE-SYN CCS for all values of RMI (y-axis) and C-FACTOR (x-axis), FIXED COUPLING 
file_names_heatmap = {['heatmap_sts_ccs_rmi_csfac_coup_' num2str(two_node_couplings(coup_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('double-syn CCS: rmi vs c-scaling (n = %d, model order = %d, coupling = %d)\n',n,p,two_node_couplings(coup_index))};
	
phiidSyn_CCS_plot = array2table(phiidSyn_CCS, 'RowNames', all_rmi_str,'VariableNames', c_factor_str);
phiidSyn_CCS_plot = {phiidSyn_CCS_plot};

plot_heatmap(phiidSyn_CCS_plot, file_names_heatmap, titles_heatmap, x_axis_csfac, ...
	y_axis_rmi, x_label_csfac, y_label_rmi, pathout_plots);
%}
