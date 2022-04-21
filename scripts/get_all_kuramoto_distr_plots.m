function get_all_kuramoto_distr_plots(data, nbins, network, A_vec, beta_vec, all_npoints, ...
		pathout_data_sim_time_series, pathout_plots_distributions) 
	
	if strcmp(data, 'raw');
		prefix = '_';
	else
		prefix = ['_' data '_'];
	end 
	
	if strcmp(data, 'raw');
		prefix2 = [];
	else
		prefix2 = [data '_'];
	end
	
	set(0,'DefaultFigureVisible','off');

	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		rng(1);
		for i = 1:length(A_vec);
			A_str = param2str(A_vec(i));
			
			for j = 1:length(beta_vec);
				beta_str = param2str(beta_vec(j));
				
				% load simulated model with given A and beta:
				% micro variables - phases, raw signal, synchronies;
				% macro variables for practical measures for causal emergence -
				% variance of synchronies & global average pairwise synchrony
				% between communities
				var1 = load([pathout_data_sim_time_series network prefix 'phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					[prefix2 'phase']);
				var2 = load([pathout_data_sim_time_series network prefix 'synchrony_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					[prefix2 'synchrony']);
				var3 = load([pathout_data_sim_time_series network prefix 'raw_signal_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					[prefix2 'raw_signal']);
				var4 = load([pathout_data_sim_time_series network prefix 'sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					[prefix2 'sigma_chi']);
				var5 = load([pathout_data_sim_time_series network prefix 'pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					[prefix2 'pair_sync']);
				var6 = load([pathout_data_sim_time_series network prefix 'mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					[prefix2 'mean_pair_sync']);
				
				% check distributions of a subset of parameters
				variables = {var1, var2, var3, var4, var5, var6};
				titles = {[prefix2 'phase'], [prefix2 'synchrony'], [prefix2 'raw signal'], [prefix2 'sigma chi'], ...
					[prefix2 'pairwise synchrony'], [prefix2 'mean pairwise synchrony']};
				filenames = {[prefix2 'phase'], [prefix2 'sync'], [prefix2 'raw_signal'], [prefix2 'sigma_chi'], ...
					[prefix2 'pair_sync'], [prefix2 'mean_pair_sync']};
				
				starting_point = npoints/2;
				end_point = (npoints/2)+100;
				
				if (((i == 1) && (j == 1)) || ((i == 3) && (j == 3)) || ((i == 7) && (j == 7)) || ((i == 10) && (j == 10)));
					
					for h = 1:length(variables);
						variable = struct2array(variables{h});
						
						num_var = size(variable, 1);
						r = randi([1 num_var],1,3);
						
						for k = 1:length(r);
							figure;
							histogram(variable(r(k),:)', nbins);
							title({['distribution of ', titles{h}]  ['var #', num2str(r(k)), ', A = ', num2str(A_vec(i)), ', beta = ', num2str(beta_vec(j))]});
							ylabel('frequency');
							xlabel(['values of ', titles{h}, ', var #', num2str(r(k))]);
							exportgraphics(gcf, [pathout_plots_distributions network '_distr_' filenames{h} '_' num2str(r(k)) '_' A_str '_' num2str(npoints) '.png']);
						end
						
						for k = 1:length(r);
							figure;
							plot(variable(r(k),starting_point:end_point), 'LineWidth', 3);
							title({['time-series of ', titles{h}]  ['var #', num2str(r(k)), ', A = ', num2str(A_vec(i)), ', beta = ', num2str(beta_vec(j))]});
							ylabel('value');
							xlabel('time');
							ylim([min(variable(r(k),starting_point:end_point))-0.3 max(variable(r(k),starting_point:end_point))+0.3]);
							xticklabels({num2str((npoints/2)) num2str(((npoints/2)+100)) num2str(((npoints/2)+200)) num2str(((npoints/2)+300)) num2str(((npoints/2)+400)) num2str(((npoints/2)+500))});
							exportgraphics(gcf, [pathout_plots_distributions network '_time_series_' filenames{h} '_' num2str(r(k)) '_' A_str '_' num2str(npoints) '.png']);
						end
						
						close all;
						
					end
				end
			end
		end
	end

end 