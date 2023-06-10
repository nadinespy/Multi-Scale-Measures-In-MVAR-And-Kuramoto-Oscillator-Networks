%% MULTI-INFO

select_measures = {'multi_info', 'phiidDoubleRed_MMI', 'phiidDoubleSyn_MMI', ...
	'phiidCE_MMI', 'phiidDC_MMI', 'phiidCD_MMI', 'phiidUC_MMI', ...
	'phiidSyn_MMI', 'phiidTransfer_MMI', 'phiidDoubleRed_CCS', ...
	'phiidDoubleSyn_CCS', 'phiidCE_CCS', 'phiidDC_CCS', 'phiidCD_CCS', ...
	'phiidUC_CCS', 'phiidSyn_CCS', 'phiidTransfer_CCS', 'DD', 'ShannonCE', ... 
	'average_corr', 'integrated_info', 'integrated_interaction', ...
	'decoder_integration', 'causal_density', 'integrated_synergy', ...
	'time_delayed_mi'};		

% 3D SCATTERPLOT ACROSS 2-NORM & RMI

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

% 2D SCATTERPLOT ACROSS 2-NORM FOR FIXED RMI

% rearrange values for one measure by having one field per measure with all
% values for different parameters
big_struct = [];
for f = 1:size(cell2mat(table2array(measures_table)),1)
	for j = 1:size(cell2mat(table2array(measures_table)),2)
		big_struct = [big_struct; cell2mat(table2array(measures_table(f,j)))];
	end
end

rmi_index = 5;
for g = 1:length(select_measures);

	BIG_array = [];
	
	for i = 1:length(all_norms);
		
		for j = rmi_index;
			
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
	title({[select_measures{g} ' multi info across 2-norm for rmi-index = ' num2str(rmi_index) ':'], 'for each pair of them, connectivities & noise correlations', 'were randmoly sampled'});
	xlabel('2-norm');
	ylabel(select_measures{g});
	ylim([0 max_measure]);
	
	location = [pathout_plots_measures, network '_' 'rmi_index' num2str(rmi_index) '_' select_measures{g}, '.png'];
	exportgraphics(gcf, location);

end

% LINE PLOT FOR AVERAGE MULTI-INFO FOR GIVEN RMI

%% DD



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
