%% PLOT ALL MEASURES

cd '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/code/analytical/scripts/measures_random_couplings'
directories = @get_mvar_directories;
directories();

n_nodes					= 8;								% number of variables in network
m_dim						= 1;								% dimension of macro variable
dim_reduction				= {'pca', 'grassmanian'};			% dim reduction method: can be 'pca' or 'grassmanian' 
time_lag_for_model		= 1;								% number of time-lags
n_samples_couplings		= 25;								% number of random samples for connectivity matrices
n_samples_noise_corrs	= 25;								% number of random samples for noise correlation matrices
seed						= 1;
network					= [num2str(n_nodes) 'mvar' ...	% model name
	'_lag' num2str(time_lag_for_model)];
spectral_radius			= 0.9;


% variation of residuals' mutual information and matrix cross_coup_mag
all_cross_coup_mag		= linspace(0.01,1,100);		% array of norm-2 values
all_cross_coup_mag_str				= {};
for t = 1:length(all_cross_coup_mag)
	all_cross_coup_mag_str{t} = num2str(all_cross_coup_mag(t));
end

all_rmi		= linspace(0.0,1,100);		% array of rmi values 
all_rmi_str	= {};
for e = 1:length(all_rmi)
	all_rmi_str{e} = num2str(all_rmi(e));
end	  

% idiosyncratic filename for this particular simulation run
filename_table = {[num2str(n_nodes) 'mvar_lag' num2str(time_lag_for_model) '_results_dd_ce_co_info_table']};
filename = {'results_all_measures'};

load([char(pathout_data_measures), char(filename_table) '.mat'],'results_table');
load([char(pathout_data_measures), network '_' char(filename) '.mat'],'results');

% full set of possible entries in table: measures are not clustered, each single measure
% and macro variable must be enumerated here
% select_measures = {'multi_info', 'phiidDoubleRed_MMI', 'phiidDoubleSyn_MMI', ...
% 	'phiidCE_MMI', 'phiidDC_MMI', 'phiidCD_MMI', 'phiidUC_MMI', ...
% 	'phiidSyn_MMI', 'phiidTransfer_MMI', 'phiidDoubleRed_CCS', ...
% 	'phiidDoubleSyn_CCS', 'phiidCE_CCS', 'phiidDC_CCS', 'phiidCD_CCS', ...
% 	'phiidUC_CCS', 'phiidSyn_CCS', 'phiidTransfer_CCS', 'DD_pca', 'ShannonCE_pca', ... 
% 	'DD_grassmanian', 'ShannonCE_grassmanian', 'co_info_pca', 'co_info_grassmanian', ...
% 	'average_corr', 'integrated_info', 'integrated_interaction', 'decoder_integration', ...
% 	'causal_density', 'integrated_synergy', 'time_delayed_mi'};

select_measures = {'DD_pca', 'ShannonCE_pca', ... 
	'DD_grassmanian', 'ShannonCE_grassmanian', 'co_info_pca', 'co_info_grassmanian'};

% select_measures = {'DD_grassmanian', 'ShannonCE_grassmanian', 'DD_pca', 'ShannonCE_pca', 'multi_info'};


%% 3D SCATTERPLOT ACROSS CROSS-COUP MAG & RMI

for g = 1:length(select_measures);
	
	BIG_array = [];
	for i = 1:length(all_cross_coup_mag);
		
		for j = 1:length(all_rmi);
			
			temp_matrix = table2array(results_table(i,j)).(select_measures{g});
			%temp_matrix = cell2mat(table2array(results_table(i,j))).(select_measures{g});
			rmi = all_rmi(j)*ones(n_samples_couplings,1);
			
			for k = 1:n_samples_noise_corrs;
				
				cross_coup_mag = table2array(results_table(i,j)).coupling_mag(:,k);
				%cross_coup_mag = cell2mat(table2array(results_table(i,j))).two_norm(:,k);
				
				big_array = [cross_coup_mag, rmi, temp_matrix(:,k)];
				BIG_array = [BIG_array; big_array];
				
			end
			
		end
		
		disp(i)
		
	end
	
	save([char(pathout_data_measures), network '_big_array_' select_measures{g} '.mat'],'results_table',  '-v7.3');

	figure;
	scatter3(BIG_array(:,1), BIG_array(:,2), BIG_array(:,3), 0.5, ".", 'blue');
	title({strrep(select_measures{g}, '_', ' '), ' across RMI and coupling magnitude: for each pair of them,', 'connectivities and noise correlations were randmoly sampled'});
	xlabel('coupling magnitude');
	ylabel('RMI');
	zlabel(strrep(select_measures{g}, '_', ' '));
	
	location = [pathout_plots_measures, network '_' select_measures{g}, '.png'];
	exportgraphics(gcf, location);

end 

%% 2D SCATTERPLOT ACROSS CROSS-COUP MAG FOR FIXED RMI

% rearrange values for one measure by having one field per measure with all
% values for different parameters
big_struct = [];
for f = 1:size(cell2mat(table2array(results_table)),1)
	for j = 1:size(cell2mat(table2array(results_table)),2)
		big_struct = [big_struct; cell2mat(table2array(results_table(f,j)))];
	end
end

rmi_index = [1 5 10];

for p = 1:length(rmi_index)
	
	for g = 1:length(select_measures);
		
		BIG_array = [];
		
		for i = 1:length(all_cross_coup_mag);
			
			for j = rmi_index(p);
				
				temp_matrix = cell2mat(table2array(results_table(i,j))).(select_measures{g});
				rmi = all_rmi(j)*ones(n_samples_couplings,1);
				
				for k = 1:n_samples_noise_corrs;
					
					cross_coup_mag = cell2mat(table2array(results_table(i,j))).two_norm(:,k);
					
					big_array = [cross_coup_mag, temp_matrix(:,k)];
					BIG_array = [BIG_array; big_array];
					
				end
				
			end
			
			%disp(i)
			
		end
		
		% extract all measure values from cells, and take maximum for ylim
		all_measure_values = [];
		for q = 1:(length(all_rmi)*length(all_cross_coup_mag));
			all_measure_values = [all_measure_values; big_struct(q,1).(select_measures{g})];
		end
		
		max_measure = max(all_measure_values(:));
		min_measure = min(all_measure_values(:));
		
		figure;
		scatter(BIG_array(:,1), BIG_array(:,2), 20, "o", 'blue', 'filled');
		title({[select_measures{g} ' across cross-coupling magnitude for rmi-index = ' num2str(rmi_index(p)) ':'], 'for each pair of them, connectivities & noise correlations', 'were randmoly sampled'});
		xlabel('cross-coupling magnitude');
		ylabel(select_measures{g});
		ylim([min_measure max_measure]);
		
		location = [pathout_plots_measures, network '_' 'rmi_index' num2str(rmi_index(p)) '_' select_measures{g}, '.png'];
		exportgraphics(gcf, location);
		
	end
end

%% LINE PLOT & HEATMAP FOR AVERAGE MEASURE FOR GIVEN RMI

%set(1,'DefaultFigureWindowStyle','docked')

	
for g = 1:length(select_measures);
	
	BIG_array = [];
	
	for i = 1:length(all_cross_coup_mag);
		
		for j = 1:length(all_rmi);
			
			mean_temp_matrix = mean(table2array(results_table(i,j)).(select_measures{g}), 'all');
			% mean_temp_matrix = mean(cell2mat(table2array(results_table(i,j))).(select_measures{g}), 'all');
			rmi = all_rmi(j);
			cross_coup_mag = mean(table2array(results_table(i,j)).coupling_mag, 'all');
			%cross_coup_mag = mean(cell2mat(table2array(results_table(i,j))).two_norm, 'all');
			
			% store values for all parameters in one big column for line plots
			big_array = [cross_coup_mag, rmi, mean_temp_matrix];
			BIG_array = [BIG_array; big_array];
			
			% store values in matrix for heatmaps
			measure(i,j) = mean_temp_matrix;
			
		end
		
		disp(i)
		
	end
	
	% LINE PLOTS: ONE LINE FOR EACH RMI
	figure;
	
	for j = 1:length(all_rmi);
		
		rmi_index = find(BIG_array(:,2)==all_rmi(j));
		
		plot3(BIG_array(rmi_index,1), BIG_array(rmi_index,2), BIG_array(rmi_index,3), 'LineWidth', 2);
		title(['average ' strrep(select_measures{g}, '_', ' ') ' across coupling magnitude for a given rmi']);
		xlabel('coupling magnitude');
		ylabel('RMI');
		zlabel(strrep(select_measures{g}, '_', ' '));
		
		hold on
	end
	
	location = [pathout_plots_measures, network '_lineplot_' select_measures{g} '.png'];
	exportgraphics(gcf, location);
	
	% HEATMAPS: CROSS-COUP MAG VS. RMI
	
	% create evenly spaced vector within range of coupling magnitudes
	min_cross_coup_mag = min(BIG_array(:,1));
	max_cross_coup_mag = max(BIG_array(:,1));
	new_cross_coup_mag_vector = round(linspace(min_cross_coup_mag,max_cross_coup_mag,length(all_rmi)),2);
	
	% get row names for table; flip order of values so that 0 starts at
	% bottom left corner in heatmap
	new_cross_coup_mag_str = {};
	for t = 1:length(new_cross_coup_mag_vector)
		new_cross_coup_mag_str{t} = num2str(new_cross_coup_mag_vector(t));
	end
	new_cross_coup_mag_str = flip(new_cross_coup_mag_str);
	
	axis_indices = [1:10:100];
	x_axis_cross_coup_mag = new_cross_coup_mag_str(axis_indices);
	y_axis_rmi = all_rmi_str(axis_indices);
	
	new_x_axis_cross_coup_mag = {x_axis_cross_coup_mag{1}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{2}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{3}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{4}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{5}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{6}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{7}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{8}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{9}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', x_axis_cross_coup_mag{10}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

	new_y_axis_rmi = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{1}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{2}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{3}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{4}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{5}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{6}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{7}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{8}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{9}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', y_axis_rmi{10}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', };

	y_label_rmi = 'RMI';
	x_label_cross_coup_mag = 'coupling magnitude';
	
	file_names_heatmap = {[network '_heatmap_coup_mag_rmi_', select_measures{g}]};
	titles_heatmap = {strrep(select_measures{g}, '_', ' ')};
	
	% flip order of rows so that lowest cross-coup mag value is in the bottom
	% left corner in heatmap and corresponding axes labels
	measure = flip(measure,1);

	new_x_axis_cross_coup_mag = flip(new_x_axis_cross_coup_mag,2);
	new_y_axis_rmi = flip(new_y_axis_rmi,2);

	%measure = array2table(measure, 'RowNames', new_cross_coup_mag_str,'VariableNames', all_rmi_str);
	measure = {measure};
	
	plot_heatmap(measure, file_names_heatmap, titles_heatmap, new_x_axis_cross_coup_mag, ...
		new_y_axis_rmi, x_label_cross_coup_mag, y_label_rmi, pathout_plots_measures);
	
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

%% CROSS-COUP MAG DISTRIBUTION

rmi_index = 100;
BIG_array = [];
for i = 1:length(all_cross_coup_mag);
	
	for j = length(all_rmi); 
		
		temp_cross_coup_mag = cell2mat(table2array(results_table(i,j))).two_norm;
		for k = 1:n_samples_noise_corrs; 
		
			BIG_array = [BIG_array; temp_cross_coup_mag(:,k)];
			
		end

	end 
	
	disp(i)
	
end 

% scatterplot
n_cross_coup_mag = 1:length(BIG_array);
figure;
scatter(n_cross_coup_mag, BIG_array, 2, ".", 'blue');
title({'cross-coupling magnitude scatterplot'});
xlabel('index'); 
ylabel('cross-coupling magnitude'); 

location = [pathout_plots_measures, network '_cross_coup_mag_scatterplot', '.png'];
exportgraphics(gcf, location);

% histogram
figure;
histogram(BIG_array, 'FaceAlpha', 0.6);
title({'cross-coupling magnitude histogram'});
xlabel('cross-coupling magnitude'); 
ylabel('frequency'); 

location = [pathout_plots_measures, network '_two_norm_histogram', '.png'];
exportgraphics(gcf, location);

close all;
