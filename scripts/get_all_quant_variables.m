function get_all_quant_variables(network, param1_vec, param2_vec, param3_vec, micro_variable_names, ...
		macro_variable_names, quantile_number, pathout_data_sim_time_series);
	
	for q = 1:length(param1_vec);
		param1_str = num2str(param1_vec(q));
		
		rng(1);
		for i = 1:length(param2_vec);
			param2_str = param2str(param2_vec(i));
			
			for j = 1:length(param3_vec)
				param3_str = param2str(param3_vec(j));
				
				fprintf('get_all_quant_variables - loop indices: npoints: %d, param1: %d, param2: %d\n', q, i, j);
			
				
				% load variables with given parameter 1 & parameter 2
				
				% load all micro variables into one struct 'micro_variables'
				for k = 1:length(micro_variable_names)
					micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' micro_variable_names{k} '_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					micro_variable_names{k}));
				end 
				
				% load all macro variables into one struct 'macro_variables'
				for k = 1:length(macro_variable_names)
					macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' macro_variable_names{k} '_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					macro_variable_names{k}));
				end 	
				
				% quantilize all micro variables, and store them in another struct 'quant_micro_variables'
				for k = 1:length(micro_variable_names)
					[quantiles_variable, quantilized_variable] = get_quantilized_variables(micro_variables(1).(micro_variable_names{k}), quantile_number);
					quant_micro_variables(1).(['quant_' micro_variable_names{k}]) = quantilized_variable;
				end
				
				% quantilize all macro variables, and store them in another struct 'quant_macro_variables'
				for k = 1:length(macro_variable_names)
					[quantiles_variable, quantilized_variable] = get_quantilized_variables(macro_variables(1).(macro_variable_names{k}), quantile_number);
					quant_macro_variables(1).(['quant_' macro_variable_names{k}]) = quantilized_variable;
				end
				
				% get all fieldnames from 'quant_macro_variables' & 'quant_macro_variables'
				quant_micro_variable_names = fieldnames(quant_micro_variables);
				quant_macro_variable_names = fieldnames(quant_macro_variables);
				
				% save all quantilized micro variables from 'quant_micro_variables' as single variables
				for k = 1:length(quant_micro_variable_names);
					% save single element from struct
					save([pathout_data_sim_time_series network '_' quant_micro_variable_names{k} '_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'-struct', 'quant_micro_variables', quant_micro_variable_names{k});
				end 
				
				% save all quantilized macro variables from 'quant_macro_variables' as single variables
				for k = 1:length(quant_macro_variable_names);
					% save single element from struct
					save([pathout_data_sim_time_series network '_' quant_macro_variable_names{k} '_' param2_str '_' param3_str '_' param1_str '.mat'], ...
					'-struct', 'quant_macro_variables', quant_macro_variable_names{k});
				end 
				
				%load([pathout_data_sim_time_series network '_' quant_micro_variable_names{k} '_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
				%	quant_micro_variable_names{k});
				
			end
		end
	end
	
end