
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_mec_mvar_rand_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In order to run this script, 
% - dirs_and_params_for_local_machine.m
% - params_phiid_mmi_rand_cross_coup_rmi.m 
% (or any other script with the prefix 'params') must be run.
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd '/media/nadinespy/NewVolume1/work/phd/projects/mec_experiments/mec_simulations/code/multivar_autoreg/rand_cross_coup_rmi/scripts'

% define model name
network = [num2str(n_nodes) 'mvar' '_lag' num2str(ar_model_order)];

% filename for this particular model & measure
filename = {[num2str(n_nodes) 'mvar_lag' num2str(ar_model_order) ...
	'_rand_cross_coup_rmi_' measure_in_filename]};
load([char(pathout_data_measures), char(filename) '.mat'], ...
	'results');

% define measure name as in table
if strcmp(measure_in_filename, 'multi_info') 
	measure_in_struct = {'MultiInfo'};
elseif strcmp(measure_in_filename, 'phiid_mmi_ce') 
	measure_in_struct = {'PhiID_CE_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_dc') 
	measure_in_struct = {'PhiID_DC_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_cd') 
	measure_in_struct = {'PhiID_CD_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_uc') 
	measure_in_struct = {'PhiID_UC_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_double_red') 
	measure_in_struct = {'PhiID_DoubleRed_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_double_syn') 
	measure_in_struct = {'PhiID_DoubleSyn_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_syn') 
	measure_in_struct = {'PhiID_Syn_MMI'};
elseif strcmp(measure_in_filename, 'phiid_mmi_transfer') 
	measure_in_struct = {'PhiID_Transfer_MMI'};
elseif strcmp(measure_in_filename, 'phiid_ccs_ce') 
	measure_in_struct = {'PhiID_CE_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_dc') 
	measure_in_struct = {'PhiID_DC_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_cd') 
	measure_in_struct = {'PhiID_CD_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_uc') 
	measure_in_struct = {'PhiID_UC_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_double_red') 
	measure_in_struct = {'PhiID_DoubleRed_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_double_syn')
	measure_in_struct = {'PhiID_DoubleSyn_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_syn') 
	measure_in_struct = {'PhiID_Syn_CCS'};
elseif strcmp(measure_in_filename, 'phiid_ccs_transfer') 
	measure_in_struct = {'PhiID_Transfer_CCS'};
elseif strcmp(measure_in_filename, 'dd_ce_co_info')
	measure_in_struct = {'DD_PCA', 'ShannonCE_PCA', 'DDGrassMin', ...
		'DDGrassMean', 'ShannonCEGrassMax', 'ShannonCEGrassMean', ...
	'CoInfoPCA', 'CoInfoGrassMin', 'CoInfoGrassMax', 'CoInfoGrassMean'};
elseif strcmp(measure_in_filename, 'integrated_info')
	measure_in_struct = {'IntegratedInfo', 'CausalDensity', ...
		'IntegratedSynergy'};
elseif strcmp(measure_in_filename, 'control')
	measure_in_struct = {'TimeDelayedMI', 'AverageCorr'};
else error('Unknown measure ''%s''', measure_in_filename);
end


%% HEATMAP FOR AVERAGE MEASURE & LINE PLOT FOR EACH RMI

% create rmi labels for plots
all_rmi_str = {};
for e = 1:length(all_rmi)
        all_rmi_str{e} = num2str(all_rmi(e));
end

% -----------------------------------------------

%% CONTINUE HERE
% -----------------------------------------------

for g = 1:length(measure_in_struct)
	
	BIG_array = [];
	
	for i = 1:length(all_cross_coup_mag)
		for j = 1:length(all_rmi)
			
			% get mean for measure across all 
			% [n_samples_cross_coup] x [n_samples_noise_corr] 
			% for given i, j
			mean_measure = mean(table2array(results_table(i,j)) ...
				.(measure_in_struct{g}), 'all');

			% get mean global coupling across all 
			% [n_samples_cross_coup] x [n_samples_noise_corr] 
			% values per i, j
			global_coup_mag = mean(table2array(results_table(i,j)) ...
				.global_coupling_mag, 'all');
			
			% store global coupling mag, rmi, and mean measure
			% in one big column (rmi column will have 
			% [length(all_rmi)] repetitions of the same value)
			big_array = [global_coup_mag, all_rmi(j), mean_measure];
			BIG_array = [BIG_array; big_array];
			
			% store mean of measure for given i, j for heatmaps
			measure(i,j) = mean_measure;
			
		end
	end
	
	% LINE PLOTS: ONE LINE FOR EACH RMI
	%{
	figure;
	
	for j = 1:length(all_rmi)
		
		rmi_index = find(BIG_array(:,2)==all_rmi(j));
		
		plot3(BIG_array(rmi_index,1), BIG_array(rmi_index,2), ...
			BIG_array(rmi_index,3), 'LineWidth', 2);
		title(['average ' strrep(measure_in_struct{g}, '_', ' ') ...
			' across coupling magnitude for a given rmi']);
		xlabel('global coupling magnitude');
		ylabel('RMI');
		zlabel(strrep(measure_in_struct{g}, '_', ' '));
		
		hold on
	end
	
	location = [pathout_plots_measures, network '_lineplot_' ...
		measure_in_struct{g} '.png'];
	exportgraphics(gcf, location);
	%}
	
	% HEATMAPS: CROSS-COUP MAG VS. RMI
	
	% create evenly spaced vector within range of global coupling 
	% magnitudes
	min_global_coup_mag = min(BIG_array(:,1));
	max_global_coup_mag = max(BIG_array(:,1));
	global_coup_mag_vector = round(linspace(min_global_coup_mag, ...
		max_global_coup_mag, length(all_cross_coup_mag)),2);
	
	% get row names for table; flip order of values so that 0 starts at
	% bottom left corner in heatmap
	global_coup_mag_str = {};
	for t = 1:length(global_coup_mag_vector)
		global_coup_mag_str{t} = num2str(global_coup_mag_vector(t));
	end
	global_coup_mag_str = flip(global_coup_mag_str);

	% flip order of rows in [measure] so that it corresponds to 
	% lowest cross-coup mag value being in the bottom
	% left corner in heatmap
	measure = flip(measure,1);
	
	% define which values on axes to display & extract them
	axis_indices = [1:10:100];
	y_axis_global_coup_mag = global_coup_mag_str(axis_indices);
	x_axis_rmi = all_rmi_str(axis_indices);
	
	% define axes ticks
	new_y_axis_global_coup_mag = {y_axis_global_coup_mag{1}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{2}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{3}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{4}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{5}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{6}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{7}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{8}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{9}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		y_axis_global_coup_mag{10}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

	new_x_axis_rmi = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{1}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{2}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{3}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{4}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{5}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{6}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{7}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{8}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{9}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ...
		x_axis_rmi{10}, ...
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', };

	% axes labels and title
	x_label_rmi = 'RMI';
	y_label_global_coup_mag = 'global coupling magnitude';
	
	filename_heatmap = {[network '_heatmap_coup_mag_rmi_', ...
		measure_in_struct{g}]};
	title_heatmap = {strrep(measure_in_struct{g}, '_', ' ')};
	
	% put [measure] in cell (needed to use plot_heatmap())
	measure = {measure};
	
	plot_heatmap(measure, ...
		filename_heatmap, ...
		title_heatmap, ...
		new_x_axis_rmi, ...
		new_y_axis_global_coup_mag, ...
		x_label_rmi, ...
		y_label_global_coup_mag, ...
		pathout_plots_measures);
	
	clear measure;
end

%% 3D SCATTERPLOT: CLUSTERS FOR DIFFERENT MEASURES

% rearrange values for one measure by having one field per measure with all
% values for different parameters
% big_column = [];
% for f = 1:length(measure_in_struct) 
% 	for k = 1:n_samples_noise_corr;
% 		one_sample = big_struct.(measure_in_struct{f})
% 		big_column = [big_column; struct2array(big_struct.(measure_in_struct{f})(:,k)];
% 	end
% end
% 
% big_column = [];
% for f = 1:length(measure_in_struct) 
% 	for k = 1:n_samples_noise_corr;
% 		one_sample = big_struct.(measure_in_struct{f})
% 		big_column = [big_column; struct2array(big_struct.(measure_in_struct{f})(:,k)];
% 	end
% end

%% CROSS-COUP MAG DISTRIBUTION

%{
rmi_index = 100;
BIG_array = [];
for i = 1:length(all_cross_coup_mag);
	
	for j = length(all_rmi); 
		
		temp_global_coup_mag = cell2mat(table2array(results_table(i,j))).two_norm;
		for k = 1:n_samples_noise_corr; 
		
			BIG_array = [BIG_array; temp_global_coup_mag(:,k)];
			
		end

	end 
	
	disp(i)
	
end 

% scatterplot
n_global_coup_mag = 1:length(BIG_array);
figure;
scatter(n_global_coup_mag, BIG_array, 2, ".", 'blue');
title({'cross-coupling magnitude scatterplot'});
xlabel('index'); 
ylabel('cross-coupling magnitude'); 

location = [pathout_plots_measures, network '_global_coup_mag_scatterplot', '.png'];
exportgraphics(gcf, location);

% histogram
figure;
histogram(BIG_array, 'FaceAlpha', 0.6);
%}

%% 3D SCATTERPLOT ACROSS CROSS-COUP MAG & RMI

%{
for g = 1:length(measure_in_struct);
	
	BIG_array = [];
	for i = 1:length(all_cross_coup_mag);
		
		for j = 1:length(all_rmi);
			
			temp_matrix = table2array(results_table(i,j)).(measure_in_struct{g});
			%temp_matrix = cell2mat(table2array(results_table(i,j))).(measure_in_struct{g});
			rmi = all_rmi(j)*ones(n_samples_global_coup,1);
			
			for k = 1:n_samples_noise_corr;
				
				global_coup_mag = table2array(results_table(i,j)).coupling_mag(:,k);
				%global_coup_mag = cell2mat(table2array(results_table(i,j))).two_norm(:,k);
				
				big_array = [global_coup_mag, rmi, temp_matrix(:,k)];
				BIG_array = [BIG_array; big_array];
				
			end
			
		end
		
		disp(i)
		
	end
	
	save([char(pathout_data_measures), network '_big_array_' measure_in_struct{g} '.mat'],'results_table',  '-v7.3');

	figure;
	scatter3(BIG_array(:,1), BIG_array(:,2), BIG_array(:,3), 0.5, ".", 'blue');
	title({strrep(measure_in_struct{g}, '_', ' '), ' across RMI and coupling magnitude: for each pair of them,', 'connectivities and noise correlations were randmoly sampled'});
	xlabel('coupling magnitude');
	ylabel('RMI');
	zlabel(strrep(measure_in_struct{g}, '_', ' '));
	
	location = [pathout_plots_measures, network '_' measure_in_struct{g}, '.png'];
	exportgraphics(gcf, location);

end 
%}

%% 2D SCATTERPLOT ACROSS CROSS-COUP MAG FOR FIXED RMI

%{
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
	
	for g = 1:length(measure_in_struct);
		
		BIG_array = [];
		
		for i = 1:length(all_cross_coup_mag);
			
			for j = rmi_index(p);
				
				temp_matrix = cell2mat(table2array(results_table(i,j))).(measure_in_struct{g});
				rmi = all_rmi(j)*ones(n_samples_global_coup,1);
				
				for k = 1:n_samples_noise_corr;
					
					global_coup_mag = cell2mat(table2array(results_table(i,j))).two_norm(:,k);
					
					big_array = [global_coup_mag, temp_matrix(:,k)];
					BIG_array = [BIG_array; big_array];
					
				end
				
			end
			
			%disp(i)
			
		end
		
		% extract all measure values from cells, and take maximum for ylim
		all_measure_values = [];
		for q = 1:(length(all_rmi)*length(all_cross_coup_mag));
			all_measure_values = [all_measure_values; big_struct(q,1).(measure_in_struct{g})];
		end
		
		max_measure = max(all_measure_values(:));
		min_measure = min(all_measure_values(:));
		
		figure;
		scatter(BIG_array(:,1), BIG_array(:,2), 20, "o", 'blue', 'filled');
		title({[measure_in_struct{g} ' across cross-coupling magnitude for rmi-index = ' num2str(rmi_index(p)) ':'], 'for each pair of them, connectivities & noise correlations', 'were randmoly sampled'});
		xlabel('cross-coupling magnitude');
		ylabel(measure_in_struct{g});
		ylim([min_measure max_measure]);
		
		location = [pathout_plots_measures, network '_' 'rmi_index' num2str(rmi_index(p)) '_' measure_in_struct{g}, '.png'];
		exportgraphics(gcf, location);
		
	end
end

%}
