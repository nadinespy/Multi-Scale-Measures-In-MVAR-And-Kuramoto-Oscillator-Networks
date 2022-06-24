function get_all_kuramoto_practCE(network, param1_vec, param2_vec, param3_vec, param4_vec, method_practCE, disc_method, ...
		micro_variable_names, macro_variable_names, pathout_data_sim_time_series, pathout_data_pract_ce)

	if ~exist('disc_method')
		disc_method = [];
	end
	
	for q = 1:length(param1_vec);
		param1_str = num2str(param1_vec(q));
		
		for z = 1:length(param2_vec);
			param2_str = num2str(param2_vec(z));
			
			rng(1);
			for i = 1:length(param3_vec);
				param3_str = param2str(param3_vec(i));
				
				for j = 1:length(param4_vec)
					param4_str = param2str(param4_vec(j));
					
					fprintf('get_all_kuramoto_practCE - loop indices: param1: %d, param2: %d, param3: %d, param4: %d\n', q, z, i, j);
					
					if strcmp(method_practCE, 'discrete');
						
						% load variables with given parameters 1, 3, & 4
				
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							micro_variables(1).([disc_method '_' micro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' disc_method '_' micro_variable_names{k} '_' param3_str '_' param4_str '_' param1_str '.mat'], ...
								[disc_method '_' micro_variable_names{k}]));
						end
						
 						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							macro_variables(1).([disc_method '_' macro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' disc_method '_' macro_variable_names{k} '_' param3_str '_' param4_str '_' param1_str '.mat'], ...
								[disc_method '_' macro_variable_names{k}]));
						end
					
					else
						
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' micro_variable_names{k} '_' param3_str '_' param4_str '_' param1_str '.mat'], ...
								micro_variable_names{k}));
						end
						
						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' macro_variable_names{k} '_' param3_str '_' param4_str '_' param1_str '.mat'], ...
								macro_variable_names{k}));
						end
						
					end 

					% ---------------------------------------------------------------------------------------------------------------------------------------
					% practical causal emergence
					% ---------------------------------------------------------------------------------------------------------------------------------------
					
					% calculate practical CE
					
					% get practical CE, DC, & CE for all combinations of micro and macro variables
					practCE = get_practCE(micro_variables, macro_variables, param2_vec(z), method_practCE);
					
					% extract practical CE, DC, CD for a given combination of micro and macro variables, 
					% and store values in separate fields in all_practCE; repeat for all combinations 
					% of variables; repeat for all values of param3_vec and param4_vec, so that fields denoting 
					% practical CE, DC, CD for each variable combination will contain 2D matrices corresponding to 
					% values of param3_vec and param4_vec. 
					
					fieldnames_practCE = fieldnames(practCE);
					
					for l = 1:length(fieldnames(practCE));
					
						all_practCE(1).(['pract_ce_' fieldnames_practCE{l,1}])(i,j) = practCE.(fieldnames_practCE{l,1}).pract_ce;
						all_practCE(1).(['pract_dc_' fieldnames_practCE{l,1}])(i,j) = practCE.(fieldnames_practCE{l,1}).pract_dc;
						all_practCE(1).(['pract_cd_' fieldnames_practCE{l,1}])(i,j) = practCE.(fieldnames_practCE{l,1}).pract_cd;
						
					end
					
					clear practCE;
					clear macro_variables;
					clear micro_variables;
					
				end
			end

			
			%% saving practical measures for different micro & macro variables
			
			% saved filenames consist of
			% network name + type of causal emergence + number of datapoints + time-lag
			save([pathout_data_pract_ce network '_all_practCE_' param1_str '_' param2_str '.mat'], ...
				'all_practCE');
			
		end
		
		clear all_practCE;
		
	end

end 