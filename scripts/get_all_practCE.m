function get_all_practCE(network, time_series_length, measure_params, model_params, ...
	micro_variable_names, macro_variable_names, ce_variable_name, pathout_data_sim_time_series, pathout_data_pract_ce)
	
	% should do separate loops for Gaussian & discrete methods

	model_param1 = model_params{1};
	model_param2 = model_params{2};
	measure_param1 = measure_params{1};
	measure_param2 = measure_params{2};
	measure_param3 = measure_params{3};	% so far unused, as Kraskov is not yet implemented for practical CE
	measure_param4 = measure_params{4};
	measure_param5 = measure_params{5};
	
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
						
						for p = 1:length(measure_param4);
							
							for o = 1:length(measure_param5);
					
								fprintf('get_all_practCE - loop indices: time_series_length: %d, measure_param1: %d, model_param1: %d, model_param2: %d, measure_param2: %d, measure_param4: %d, measure_param5: %d\n', q, z, i, j, b, p, o);
								
								if strcmp(lower(measure_param2(b)), 'discrete');
									
									% load variables with given parameters
									
									% load all micro variables into one struct 'micro_variables'
									for k = 1:length(micro_variable_names)
										micro_variables(1).([measure_param4{p} num2str(measure_param5(o)) '_' micro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' measure_param4{p} num2str(measure_param5(o)) '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
											[measure_param4{p} num2str(measure_param5(o)) '_' micro_variable_names{k}]));
									end
									
									% load all macro variables into one struct 'macro_variables'
									for k = 1:length(macro_variable_names)
										macro_variables(1).([measure_param4{p} num2str(measure_param5(o)) '_' macro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' measure_param4{p} num2str(measure_param5(o)) '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
											[measure_param4{p} num2str(measure_param5(o)) '_' macro_variable_names{k}]));
									end
									
								else
									
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
									
								end
								
								% ---------------------------------------------------------------------------------------------------------------------------------------
								% practical causal emergence
								% ---------------------------------------------------------------------------------------------------------------------------------------
								
								% calculate practical CE
								
								% get practical CE, DC, & CE for all combinations of micro and macro variables
								practCE = get_practCE(micro_variables, macro_variables, measure_param1(z), measure_param2{b});
								
								% extract practical CE, DC, CD for a given combination of micro and macro variables,
								% and store values in separate fields in all_practCE; repeat for all combinations
								% of variables; repeat for all values of model_param1 and model_param2, so that fields denoting
								% practical CE, DC, CD for each variable combination will contain 2D matrices corresponding to
								% values of model_param1 and model_param2.
								
								fieldnames_practCE = fieldnames(practCE);

								for l = 1:length(fieldnames(practCE));
						
									temp_practCE(1).(['pract_ce_' fieldnames_practCE{l,1}])(i,j) = practCE.(fieldnames_practCE{l,1}).pract_ce;
									temp_practCE(1).(['pract_dc_' fieldnames_practCE{l,1}])(i,j) = practCE.(fieldnames_practCE{l,1}).pract_dc;
									temp_practCE(1).(['pract_cd_' fieldnames_practCE{l,1}])(i,j) = practCE.(fieldnames_practCE{l,1}).pract_cd;
									
								end
								
								fieldnames_temp_practCE = fieldnames(temp_practCE);
								
								if strcmp(lower(measure_param2(b)), 'discrete');	
									% fieldnames consist of measure (practical CE, DC, or CD) + micro variable (with bin number) + macro variable (with bin number)
									for l = 1:length(fieldnames_temp_practCE);
										new_fieldnames_temp_practCE{l} = [fieldnames_temp_practCE{l}];
									end
								else 
									% fieldnames consist of measure (practical CE, DC, or CD) + micro variable + macro variable + pract CE method
									for l = 1:length(fieldnames_temp_practCE);
										new_fieldnames_temp_practCE{l} = [fieldnames_temp_practCE{l} '_gauss'];
									end
								end
								
								for l = 1:length(new_fieldnames_temp_practCE);
									all_practCE(1).(new_fieldnames_temp_practCE{l}) = temp_practCE(1).(fieldnames_temp_practCE{l});
								end
								
								clear practCE;
								clear macro_variables;
								clear micro_variables;
								
							end
						end		
					end
				end
			end

			
			%% saving practical measures for different micro & macro variables
			
			% saved filenames consist of
			% network name + type of causal emergence + number of datapoints + time-lag
			save([pathout_data_pract_ce network '_all_practCE_' ce_variable_name '_' time_series_length_str '_' measure_param1_str '.mat'], ...
				'all_practCE');
			
		end
		
		clear all_practCE;
		
	end
	
end 
