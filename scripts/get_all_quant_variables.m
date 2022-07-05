function get_all_quant_variables(network, time_series_length, measure_params, model_params, micro_variable_names, ...
		macro_variable_names, pathout_data_sim_time_series);
	
	model_param1 = model_params{1};
	model_param2 = model_params{2};
	measure_param5 = measure_params{5};
	
	for q = 1:length(time_series_length);
		time_series_length_str = num2str(time_series_length(q));
		
		rng(1);
		for i = 1:length(model_param1);
			model_param1_str = param2str(model_param1(i));
			
			for j = 1:length(model_param2)
				model_param2_str = param2str(model_param2(j));
				
				for o = 1:length(measure_param5);
				
					fprintf('get_all_quant_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d, measure_param5: %d\n', q, i, j, o);
					
					
					% load variables with given parameter 1 & parameter 2
					
					% load all micro variables into one struct 'micro_variables'
					for k = 1:length(micro_variable_names)
						micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
							micro_variable_names{k}));
					end
					
					% load all macro variables into one struct 'macro_variables'
					for k = 1:length(macro_variable_names)
						macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
							macro_variable_names{k}));
					end
					
					% quantilize all micro variables, and store them in another struct 'quant_micro_variables'
					for k = 1:length(micro_variable_names)
						[quantiles_variable, quantilized_variable] = get_quantilized_variables(micro_variables(1).(micro_variable_names{k}), measure_param5(o));
						quant_micro_variables(1).(['quant_' micro_variable_names{k}]) = quantilized_variable;
					end
					
					% quantilize all macro variables, and store them in another struct 'quant_macro_variables'
					for k = 1:length(macro_variable_names)
						[quantiles_variable, quantilized_variable] = get_quantilized_variables(macro_variables(1).(macro_variable_names{k}), measure_param5(o));
						quant_macro_variables(1).(['quant_' macro_variable_names{k}]) = quantilized_variable;
					end
					
					% get all fieldnames from 'quant_macro_variables' & 'quant_macro_variables'
					quant_micro_variable_names = fieldnames(quant_micro_variables);
					quant_macro_variable_names = fieldnames(quant_macro_variables);
					
					% save all quantilized micro variables from 'quant_micro_variables' as single variables
					% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints + number of quantiles
					for k = 1:length(quant_micro_variable_names);
						% save single element from struct
						save([pathout_data_sim_time_series network '_' quant_micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '_' num2str(measure_param5(o)) '.mat'], ...
							'-struct', 'quant_micro_variables', quant_micro_variable_names{k});
					end
					
					% save all quantilized macro variables from 'quant_macro_variables' as single variables
					% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints + number of quantiles
					for k = 1:length(quant_macro_variable_names);
						% save single element from struct
						save([pathout_data_sim_time_series network '_' quant_macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '_' num2str(measure_param5(o)) '.mat'], ...
							'-struct', 'quant_macro_variables', quant_macro_variable_names{k});
					end
					
					%load([pathout_data_sim_time_series network '_' quant_micro_variable_names{k} '_' A_str '_' beta_str  '_' num2str(npoints) '.mat'], ...
					%	quant_micro_variable_names{k});
					
				end
			end
		end
	end
	
end
