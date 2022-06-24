function get_all_disc_variables(network, param1_vec, param2_vec, all_npoints, micro_variable_names, ...
		macro_variable_names, bin_number, pathout_data_sim_time_series);
	
	for q = 1:length(all_npoints);
		npoints = all_npoints(q);
		
		rng(1);
		for i = 1:length(param1_vec);
			param1_str = param2str(param1_vec(i));
			
			for j = 1:length(param2_vec)
				param2_str = param2str(param2_vec(j));
				
				fprintf('get_all_disc_variables - loop indices: npoints: %d, param1: %d, param2: %d\n', q, i, j);
			
				
				% load variables with given parameter 1 & parameter 2
				
				% load all micro variables into one struct 'micro_variables'
				for k = 1:length(micro_variable_names)
					micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' micro_variable_names{k} '_' param1_str '_' param2_str '_' num2str(npoints) '.mat'], ...
					micro_variable_names{k}));
				end 
				
				% load all macro variables into one struct 'macro_variables'
				for k = 1:length(macro_variable_names)
					macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' macro_variable_names{k} '_' param1_str '_' param2_str '_' num2str(npoints) '.mat'], ...
					macro_variable_names{k}));
				end 	
				
				% discretize all micro variables, and store them in another struct 'disc_micro_variables'
				for k = 1:length(micro_variable_names)
					discretized_variable = get_discretized_variables(micro_variables(1).(micro_variable_names{k}), bin_number);
					disc_micro_variables(1).(['disc_' micro_variable_names{k}]) = discretized_variable;
				end
				
				% discretize all macro variables, and store them in another struct 'disc_macro_variables'
				for k = 1:length(macro_variable_names)
					discretized_variable = get_discretized_variables(macro_variables(1).(macro_variable_names{k}), bin_number);
					disc_macro_variables(1).(['disc_' macro_variable_names{k}]) = discretized_variable;
				end
				
				% get all fieldnames from 'disc_macro_variables' & 'disc_macro_variables'
				disc_micro_variable_names = fieldnames(disc_micro_variables);
				disc_macro_variable_names = fieldnames(disc_macro_variables);
				
				% save all discretized micro variables from 'disc_micro_variables' as single variables
				for k = 1:length(disc_micro_variable_names);
					% save single element from struct
					save([pathout_data_sim_time_series network '_' disc_micro_variable_names{k} '_' param1_str '_' param2_str  '_' num2str(npoints) '.mat'], ...
					'-struct', 'disc_micro_variables', disc_micro_variable_names{k});
				end 
				
				% save all discretized macro variables from 'disc_macro_variables' as single variables
				for k = 1:length(disc_macro_variable_names);
					% save single element from struct
					save([pathout_data_sim_time_series network '_' disc_macro_variable_names{k} '_' param1_str '_' param2_str  '_' num2str(npoints) '.mat'], ...
					'-struct', 'disc_macro_variables', disc_macro_variable_names{k});
				end 
				
				%load([pathout_data_sim_time_series network '_' disc_micro_variable_names{k} '_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
				%	disc_micro_variable_names{k});
				
			end
		end
	end
	
end