%% PLOTTING: INTEGRATION IN TWO-NODE MVAR NETWORKS - NOISE CORRELATION & COUPLING, FIXED C-FACTOR

%% line plots

% {
% plot INTEGRATION (y-axis) and COUPLING (x-axis), one line per NOISE CORRELATION, FIXED C-FACTOR
figure(1); clf
plot(noise_corrs',I');

lgd = legend(num2str(noise_corrs'));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs coupling (n = %d, model order = %d, c-factor = %d),\n single line for each noise correlation', ...
	n,time_lag,csfac(csfac_index)));
xlabel('coupling');
ylabel('integration (nats)');

grid on

file_names = {['line_noise_corr_coup_csfac_' num2str(csfac(csfac_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);

% plot INTEGRATION (y-axis) and NOISE CORRELATION (x-axis), one line per COUPLINGS, FIXED C-FACTOR
figure(1); clf
plot(noise_corrs',I);

lgd = legend(num2str(two_node_couplings'));
lgd.FontSize = 5;
legend('Location','eastoutside');
title(sprintf('Integration vs noise correlation (n = %d, model order = %d, c-factor = %d),\n single line for each coupling', ...
	n,time_lag,csfac(csfac_index)));
xlabel('noise correlation');
ylabel('integration (nats)');

grid on

file_names = {['line_coup_noise_corr_csfac_' num2str(csfac(csfac_index)) '_' signs_noise_corrs(1:3)]};
file_names = strrep(file_names,'.','');
location = string(strcat(pathout_plots, file_names, '.png'));
exportgraphics(gcf, location);
%}

%% heatmaps

% {
% -------------------------------------------------------------------------
% plot INTEGRATION for all values of NOISE CORRELATIONS (y-axis) and COUPLINGS (x-axis), FIXED C-FACTOR
y_axis_noise_corr = {'0.0', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.32', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.56', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.9'};

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
noise_corr_str = {};
for t = 1:length(noise_corrs)
	noise_corr_str{t} = num2str(noise_corrs(t));
end

y_label_noise_corr = 'noise correlation';
x_label_coup = 'coupling';

file_names_heatmap = {['heatmap_noise_corr_coup_csfac_' num2str(csfac(csfac_index)) '_' signs_noise_corrs(1:3)]};
file_names_heatmap = strrep(file_names_heatmap,'.','');
titles_heatmap = {sprintf('noise correlation vs coupling (n = %d, model order = %d, c-factor = %d)\n', ...
	n,time_lag,csfac(csfac_index))};
	
integration = array2table(I, 'RowNames', noise_corr_str,'VariableNames', two_node_coupling_str);
integration = {integration};

plot_heatmap(integration, file_names_heatmap, titles_heatmap, x_axis_coup, ...
	y_axis_noise_corr, x_label_coup, y_label_noise_corr, pathout_plots);
%}