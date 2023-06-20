%% PLOT ALL MEASURES

load([char(pathout_data_measures), network '_measures_table.mat'],'measures_table');
load([char(pathout_data_measures), network '_measures.mat'],'measures');

select_measures = {'multi_info', 'phiidDoubleRed_MMI', 'phiidDoubleSyn_MMI', ...
	'phiidCE_MMI', 'phiidDC_MMI', 'phiidCD_MMI', 'phiidUC_MMI', ...
	'phiidSyn_MMI', 'phiidTransfer_MMI', 'phiidDoubleRed_CCS', ...
	'phiidDoubleSyn_CCS', 'phiidCE_CCS', 'phiidDC_CCS', 'phiidCD_CCS', ...
	'phiidUC_CCS', 'phiidSyn_CCS', 'phiidTransfer_CCS', 'DD_pca', 'ShannonCE_pca', ... 
	'DD_grassmanian', 'ShannonCE_grassmanian', 'average_corr', 'integrated_info', ...
	'integrated_interaction', 'decoder_integration', 'causal_density', ...
	'integrated_synergy', 'time_delayed_mi'};		

%% 3D SCATTERPLOT ACROSS 2-NORM & RMI

for g = 1:length(select_measures);
	
	BIG_array = [];
	for i = 1:length(all_norms);
		
		for j = 1:length(all_rmi);
			
			temp_matrix = cell2mat(table2array(measures_table(i,j))).(select_measures{g});
			rmi = all_rmi(j)*ones(n_samples_couplings,1);
			
			for k = 1:n_samples_noise_corrs;
				
				norm = cell2mat(table2array(measures_table(i,j))).two_norm(:,k);
				
				big_array = [norm, rmi, temp_matrix(:,k)];
				BIG_array = [BIG_array; big_array];
				
			end
			
		end
		
		%disp(i)
		
	end
	
	figure;
	scatter3(BIG_array(:,1), BIG_array(:,2), BIG_array(:,3), 10, ".", 'blue');
	title({select_measures{g} ' across RMI and 2-norm: for each pair of them,', 'connectivities and noise correlations were randmoly sampled'});
	xlabel('2-norm');
	ylabel('RMI');
	zlabel(select_measures{g});
	
	location = [pathout_plots_measures, network '_' select_measures{g}, '.png'];
	exportgraphics(gcf, location);

end 

%% 2D SCATTERPLOT ACROSS 2-NORM FOR FIXED RMI

% rearrange values for one measure by having one field per measure with all
% values for different parameters
big_struct = [];
for f = 1:size(cell2mat(table2array(measures_table)),1)
	for j = 1:size(cell2mat(table2array(measures_table)),2)
		big_struct = [big_struct; cell2mat(table2array(measures_table(f,j)))];
	end
end

rmi_index = [1 5 10];

for p = 1:length(rmi_index)
	
	for g = 1:length(select_measures);
		
		BIG_array = [];
		
		for i = 1:length(all_norms);
			
			for j = rmi_index(p);
				
				temp_matrix = cell2mat(table2array(measures_table(i,j))).(select_measures{g});
				rmi = all_rmi(j)*ones(n_samples_couplings,1);
				
				for k = 1:n_samples_noise_corrs;
					
					norm = cell2mat(table2array(measures_table(i,j))).two_norm(:,k);
					
					big_array = [norm, temp_matrix(:,k)];
					BIG_array = [BIG_array; big_array];
					
				end
				
			end
			
			%disp(i)
			
		end
		
		% extract all measure values from cells, and take maximum for ylim
		all_measure_values = [];
		for q = 1:(length(all_rmi)*length(all_norms));
			all_measure_values = [all_measure_values; big_struct(q,1).(select_measures{g})];
		end
		
		max_measure = max(all_measure_values(:));
		
		figure;
		scatter(BIG_array(:,1), BIG_array(:,2), 20, "o", 'blue', 'filled');
		title({[select_measures{g} ' across 2-norm for rmi-index = ' num2str(rmi_index(p)) ':'], 'for each pair of them, connectivities & noise correlations', 'were randmoly sampled'});
		xlabel('2-norm');
		ylabel(select_measures{g});
		ylim([0 max_measure]);
		
		location = [pathout_plots_measures, network '_' 'rmi_index' num2str(rmi_index(p)) '_' select_measures{g}, '.png'];
		exportgraphics(gcf, location);
		
	end
end

%% LINE PLOT & HEATMAP FOR AVERAGE MEASURE FOR GIVEN RMI

set(1,'DefaultFigureWindowStyle','docked')

	
for g = 1:length(select_measures);
	
	BIG_array = [];
	
	for i = 1:length(all_norms);
		
		for j = 1:length(all_rmi);
			
			mean_temp_matrix = mean(cell2mat(table2array(measures_table(i,j))).(select_measures{g}), 'all');
			rmi = all_rmi(j);
			norm = mean(cell2mat(table2array(measures_table(i,j))).two_norm, 'all');
			
			% store values for all parameters in one big column for line plots
			big_array = [norm, rmi, mean_temp_matrix];
			BIG_array = [BIG_array; big_array];
			
			% store values in matrix for heatmaps
			measure(i,j) = mean_temp_matrix;
			
		end
		
		%disp(i)
		
	end
	
	% LINE PLOTS: ONE LINE FOR EACH RMI
	figure;
	
	for j = 1:length(all_rmi);
		
		rmi_index = find(BIG_array(:,2)==all_rmi(j));
		
		plot3(BIG_array(rmi_index,1), BIG_array(rmi_index,2), BIG_array(rmi_index,3));
		title(['average ' select_measures{g} ' across 2-norm for a given rmi']);
		xlabel('2-norm');
		ylabel('RMI');
		zlabel(select_measures{g});
		
		hold on
	end
	
	location = [pathout_plots_measures, network '_lineplot_' select_measures{g} '.png'];
	exportgraphics(gcf, location);
	
	% HEATMAPS: 2-NORM VS. RMI
	
	% create evenly spaced vector within range of 2-norms
	min_norm = min(BIG_array(:,1));
	max_norm = max(BIG_array(:,1));
	new_norm_vector = round(linspace(min_norm,max_norm,length(all_rmi)),2);
	
	% get row names for table; flip order of values so that 0 starts at
	% bottom left corner in heatmap
	new_norm_str = {};
	for t = 1:length(new_norm_vector)
		new_norm_str{t} = num2str(new_norm_vector(t));
	end
	new_norm_str = flip(new_norm_str);
	
	y_axis_rmi = new_norm_str;
	x_axis_norm = all_rmi_str;
	
	y_label_rmi = '2-norm';
	x_label_norm = 'RMI';
	
	file_names_heatmap = {[network '_heatmap_norm_rmi_', select_measures{g}]};
	titles_heatmap = {select_measures{g}};
	
	% flip order of rows so that lowest 2-norm value is in the bottom
	% left corner in heatmap
	measure = flip(measure,1);
	measure = array2table(measure, 'RowNames', new_norm_str,'VariableNames', all_rmi_str);
	measure = {measure};
	
	plot_heatmap(measure, file_names_heatmap, titles_heatmap, x_axis_norm, ...
		y_axis_rmi, x_label_norm, y_label_rmi, pathout_plots_measures);
	
	clear measure;
end



%% 3D SCATTERPLOT: CLUSTERS FOR DIFFERENT MEASURES

% rearrange values for one measure by having one field per measure with all
% values for different parameters
% big_column = [];
% for f = 1:length(select_measures) 
% 	for k = 1:n_samples_noise_corrs;
% 		one_sample = big_struct.(select_measures{f})
% 		big_column = [big_column; struct2array(big_struct.(select_measures{f})(:,k)];
% 	end
% end
% 
% big_column = [];
% for f = 1:length(select_measures) 
% 	for k = 1:n_samples_noise_corrs;
% 		one_sample = big_struct.(select_measures{f})
% 		big_column = [big_column; struct2array(big_struct.(select_measures{f})(:,k)];
% 	end
% end

%% 2-NORM DISTRIBUTION

rmi_index = 100;
BIG_array = [];
for i = 1:length(all_norms);
	
	for j = length(all_rmi); 
		
		temp_norm = cell2mat(table2array(measures_table(i,j))).two_norm;
		for k = 1:n_samples_noise_corrs; 
		
			BIG_array = [BIG_array; temp_norm(:,k)];
			
		end

	end 
	
	disp(i)
	
end 

% scatterplot
n_norms = 1:length(BIG_array);
figure;
scatter(n_norms, BIG_array, 2, ".", 'blue');
title({'2-norm scatterplot'});
xlabel('index'); 
ylabel('2-norm'); 

location = [pathout_plots_measures, network '_two_norm_scatterplot', '.png'];
exportgraphics(gcf, location);

% histogram
figure;
histogram(BIG_array, 'FaceAlpha', 0.6);
title({'2-norm histogram'});
xlabel('2-norm'); 
ylabel('frequency'); 

location = [pathout_plots_measures, network '_two_norm_histogram', '.png'];
exportgraphics(gcf, location);

close all;
