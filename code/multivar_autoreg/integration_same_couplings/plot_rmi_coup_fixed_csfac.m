%% PLOTTING: INTEGRATION IN TWO-NODE MVAR NETWORKS - RMI & COUPLING, FIXED C-FACTOR

%% line plots

% {
% plot INTEGRATION (y-axis) and COUPLING (x-axis), one line per RMI, FIXED C-FACTOR
figure(1); clf
plot(two_node_couplings',I');

lgd = legend(num2str(all_rmi'));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs coupling (n = %d, model order = %d, c-factor = %d),\n single line for each rmi', ...
	n,time_lag,csfac(csfac_index)));
xlabel('coupling');
ylabel('integration (nats)');

grid on

file_names = {['line_rmi_coup_csfac_' num2str(csfac(csfac_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);

% plot INTEGRATION (y-axis) and RMI (x-axis), one line per COUPLING, FIXED C-FACTOR
figure(1); clf
plot(all_rmi',I);

lgd = legend(num2str(two_node_couplings'));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs rmi (n = %d, model order = %d, c-factor = %d),\n single line for each coupling', ...
	n,time_lag,csfac(csfac_index)));
xlabel('rmi');
ylabel('integration (nats)');

grid on

file_names = {['line_coup_rmi_csfac_' num2str(csfac(csfac_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);
%}

%% heatmaps

% {
% -------------------------------------------------------------------------
% plot INTEGRATION for all values of RMI (y-axis) and COUPLING (x-axis), FIXED C-FACTOR 
y_axis_rmi = {'0.0', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.33', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.66', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '1'};

x_axis_coup = {'0.0045', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.15', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.45'};

% get column names for table
two_node_coupling_str = {};
for t = 1:length(two_node_couplings)
	two_node_coupling_str{t} = num2str(two_node_couplings(t));
end

% get row names for table
all_rmi_str = {};
for e = 1:length(all_rmi)
	all_rmi_str{e} = num2str(all_rmi(e));
end

y_label_rmi	= 'rmi';
x_label_coup= 'coupling';

file_name_heatmap = {['heatmap_rmi_coup_csfac_' num2str(csfac(csfac_index)) '_' signs_noise_corrs(1:3)]};
file_name_heatmap = strrep(file_name_heatmap,'.','');
title_heatmap = {sprintf('rmi vs coupling (n = %d, model order = %d, c-factor = %d)\n', ...
	n,time_lag,csfac(csfac_index))};
	
integration = array2table(I, 'RowNames', all_rmi_str,'VariableNames', two_node_coupling_str);
integration = {integration};

plot_heatmap(integration, file_name_heatmap, title_heatmap, x_axis_coup, ...
	y_axis_rmi, x_label_coup, y_label_rmi, pathout_plots);
%}