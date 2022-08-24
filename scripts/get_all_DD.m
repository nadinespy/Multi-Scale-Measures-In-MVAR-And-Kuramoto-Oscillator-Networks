function get_all_DD(network, time_series_length, measure_params, measure_params_dd, model_params, micro_variable_names, ...
	macro_variable_names, dd_variable_name, pathout_data_sim_time_series, pathout_data_dd)

	% not a good solution, as order of variables is hard-coded - doesn't allow for 'measure_params' to *not* contain, 
	% e. g., kraskov parameter
	model_param1 = model_params{1};
	model_param2 = model_params{2};
	measure_param1 = measure_params{1};
	measure_param2 = measure_params{2};
	measure_param3 = measure_params{3};	
	measure_param4 = measure_params{4};
	measure_param5 = measure_params{5};
	
	measure_param1_dd = measure_params_dd{1};
		
	for q = 1:length(time_series_length);
		time_series_length_str = num2str(time_series_length(q));
		
		for z = 1:length(measure_param1);
			measure_param1_str = num2str(measure_param1(z));
			
			rng(1);
			for i = 1:length(model_param1);
				model_param1_str = param2str(model_param1(i));
				
				for j = 1:length(model_param2)
					model_param2_str = param2str(model_param2(j));
					
					for b = 1:length(measure_param2);
						
						if strcmp(lower(measure_param2{b}), 'kraskov');
							
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
							
							for c = 1:length(measure_param3)
								
								for d = 1:length(measure_param1_dd)
									
									if measure_param1_dd(d) <= measure_param1(z)

										% ---------------------------------------------------------------------------------------------------------------------------------------
										% dynamical dependence
										% ---------------------------------------------------------------------------------------------------------------------------------------
										
										% calculate dynamical independence (DD)
										
										fprintf('get_all_DD - loop indices: time_series_length: %d, measure_param1: %d, model_param1: %d, model_param2: %d, measure_param2: %d, measure_param3: %d, measure_param1_dd: %d\n', q, z, i, j, b, c, d);
										
										% get DD for all combinations of micro and top-level macro variables
										DD = get_DD(micro_variables, macro_variables, measure_param2{b}, measure_param1(z), measure_param1_dd(d), measure_param3(c));
										
										fieldnames_DD = fieldnames(DD);
										
										for l = 1:length(fieldnames(DD));
											temp_DD(1).(['dd_' fieldnames_DD{l,1}])(i,j) = DD.(fieldnames_DD{l,1});
										end
										
										fieldnames_temp_DD = fieldnames(temp_DD);
										
										% fieldnames consist of measure (practical CE, DC, or CD) + micro variable + macro variable + kraskov number + time-step
										for l = 1:length(fieldnames_temp_DD);
											new_fieldnames_temp_DD{l} = [fieldnames_temp_DD{l} '_krask' num2str(measure_param3(c)) '_tstep' num2str(measure_param1_dd(d))];
										end
										
										for l = 1:length(new_fieldnames_temp_DD);
											all_DD(1).(new_fieldnames_temp_DD{l}) = temp_DD(1).(fieldnames_temp_DD{l});
										end
									else 
										continue 
									end
								end 
							end 
							
							clear DD;
							clear macro_variables;
							clear micro_variables;
							
						elseif strcmp(lower(measure_param2{b}), 'discrete') 
						
							for p = 1:length(measure_param4);
								
								for o = 1:length(measure_param5);
									
									% load all micro variables into one struct 'micro_variables'
									for k = 1:length(micro_variable_names)
										micro_variables(1).([measure_param4{p} '_' micro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' measure_param4{p} num2str(measure_param5(o)) '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
											[measure_param4{p} num2str(measure_param5(o)) '_' micro_variable_names{k}]));
									end
									
									% load all macro variables into one struct 'macro_variables'
									for k = 1:length(macro_variable_names)
										macro_variables(1).([measure_param4{p} '_' macro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' measure_param4{p} num2str(measure_param5(o)) '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
											[measure_param4{p} num2str(measure_param5(o)) '_' macro_variable_names{k}]));
									end
									
									for u = 1:length(measure_param1_dd)
					
										fprintf('get_all_DD - loop indices: time_series_length: %d, measure_param1: %d, model_param1: %d, model_param2: %d, measure_param2: %d, measure_param4: %d, measure_param5: %d, measure_param1_dd: %d\n', q, z, i, j, b, p, o, u);
										
										% ---------------------------------------------------------------------------------------------------------------------------------------
										% dynamical dependence
										% ---------------------------------------------------------------------------------------------------------------------------------------
										
										% calculate dynamical independence (DD)
										
										% get DD for all combinations of micro and top-level macro variables
										
										DD = get_DD(micro_variables, macro_variables, measure_param2{b}, measure_param1(z), measure_param1_dd(u));
										
										fieldnames_DD = fieldnames(DD);
										
										for l = 1:length(fieldnames(DD));
											temp_DD(1).(['dd_' fieldnames_DD{l,1}])(i,j) = DD.(fieldnames_DD{l,1});
										end
										
										fieldnames_temp_DD = fieldnames(temp_DD);
										
										% fieldnames consist of measure (practical CE, DC, or CD) + micro variable + macro variable + disc number
 
										for l = 1:length(fieldnames_temp_DD);
											new_fieldnames_temp_DD{l} = [fieldnames_temp_DD{l} '_tstep' num2str(measure_param1_dd(d))];
										end
									
										for l = 1:length(new_fieldnames_temp_DD);
											all_DD(1).(new_fieldnames_temp_DD{l}) = temp_DD(1).(fieldnames_temp_DD{l});
										end
										
										clear DD;
										clear macro_variables;
										clear micro_variables;
										
									end
								end
							end
							
						elseif strcmp(lower(measure_param2{b}), 'gaussian')
							
							for u = 1:length(measure_param1_dd)
								
								fprintf('get_all_DD - loop indices: time_series_length: %d, measure_param1: %d, model_param1: %d, model_param2: %d, measure_param2: %d, measure_param4: %d, measure_param5: %d, measure_param1_dd: %d\n', q, z, i, j, b, p, u);	

								% get DD for all combinations of micro and top-level macro variables
								DD = get_DD(micro_variables, macro_variables, measure_param2{b}, measure_param1(z), measure_param1_dd);
								
								fieldnames_DD = fieldnames(DD);
								
								for l = 1:length(fieldnames(DD));
									temp_DD(1).(['dd_' fieldnames_DD{l,1}])(i,j) = DD.(fieldnames_DD{l,1});
								end
								
								fieldnames_temp_DD = fieldnames(temp_DD);
								
								% fieldnames consist of measure (practical CE, DC, or CD) + micro variable + macro variable + kraskov number + time-step
								for l = 1:length(fieldnames_temp_DD);
									new_fieldnames_temp_DD{l} = [fieldnames_temp_DD{l} '_gauss_tstep' num2str(measure_param1_dd(d))];
								end
								
								for l = 1:length(new_fieldnames_temp_DD);
									all_DD(1).(new_fieldnames_temp_DD{l}) = temp_DD(1).(fieldnames_temp_DD{l});
								end
							end
							
							clear DD;
							clear macro_variables;
							clear micro_variables;

						end
					end
				end
			end

			
			%% saving practical measures for different micro & macro variables
			
			% saved filenames consist of
			% network name + type of causal emergence + number of datapoints + time-lag
			save([pathout_data_dd network '_all_DD_' dd_variable_name '_' num2str(time_series_length(q)) '_' num2str(measure_param1(z)) '.mat'], ...
				'all_DD');
			
		end
		
		clear all_DD;
		
	end

end 
								
							
