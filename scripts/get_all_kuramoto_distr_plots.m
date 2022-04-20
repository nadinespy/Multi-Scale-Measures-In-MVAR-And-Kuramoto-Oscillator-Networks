function get_all_kuramoto_distr_plots(nbins, network, A_vec, beta_vec, all_npoints, ...
		pathout_data_sim_time_series, pathout_plots_distributions) 
	
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
				load([pathout_data_sim_time_series network '_phase_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'phase');
				load([pathout_data_sim_time_series network '_synchrony_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'synchrony');
				load([pathout_data_sim_time_series network '_raw_signal_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'raw_signal');
				load([pathout_data_sim_time_series network '_sigma_chi_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'sigma_chi');
				load([pathout_data_sim_time_series network '_mean_pair_sync_' A_str '_' beta_str '_' num2str(npoints) '.mat'], ...
					'grand_mean_pair_sync');
				
				% check distributions of a subset of parameters
				variables = {phase, raw_signal, synchrony, sigma_chi, grand_mean_pair_sync};
				titles = {'phase', 'raw signal', 'synchrony', 'sigma chi', 'global mean pairwise synchrony'};
				filenames = {'phase', 'raw_signal', 'sync', 'sigma_chi', 'pair_sync'};
				
				starting_point = npoints/2;
				end_point = (npoints/2)+100;
				
				if (((i == 1) && (j == 1)) || ((i == 3) && (j == 3)) || ((i == 7) && (j == 7)) || ((i == 10) && (j == 10)));
					
					for h = 1:length(variables);
						variable = variables{h};
						
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