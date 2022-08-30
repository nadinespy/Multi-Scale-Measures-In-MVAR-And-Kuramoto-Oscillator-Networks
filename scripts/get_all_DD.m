function output_struct = get_all_DD(input_struct, ...
		model_params1, ...
		model_params2, ...
		measure, ...
		time_length, ...
		method, ...
		time_lag, ...
		time_steps, ...
		micro_variable_names, ...
		macro_variable_names, ...
		network, ...
		pathout_data_sim_time_series, ...
		kraskov_param, ...
		disc_method, ...
		bin_number);
	
	% Inputs:	input_struct			struct with fields 
	%
	%							'time_length',
	%							'measure', 
	%							'method', 
	%							'time_lag', 
	%							'disc_method', 
	%							'bin_number', 
	%							'kraskov_param',
	%							'time_step',
	%							'red_func',
	%							'results'
	%
	%							where 'results' contains
	%							fieldnames according to
	%							micro-macro combinations,
	%							each of which contain a 
	%							table with model_params1
	%							and model_params2 as rows/
	%							columns, initialized with 
	%							zeros; all other fields
	%							are empty
	% 
	%		model_params1			float array
	%		model_params2			float array
	%		measure				string
	%		time_length				int
	%		method				string
	%		time_lag				int
	%		time_steps				int array
	%		micro_variable_names		cell array with chars
	%		macro_variable_names		cell array with chars
	%		network				string
	%		pathout_data_sim_time_series	string
	%		kraskov_param			int
	%		disc_method				string
	%		bin_number				int	
	%
	% Outputs: output_struct			same struct as input_struct,
	%							but with DD values (as 
	%							opposed to zeros) for 
	%							different micro/macro 
	%							combinations, and different
	%							model parameters in 
	%							output_struct.results
	
	output_struct = input_struct;
	fieldnames_results = fieldnames(input_struct(1,1).results(1,1));
	time_length_str = num2str(time_length);
	
	% loop over time_steps, model_params1, model_params2
	for a = 1:length(time_steps)
		time_step = time_steps(a);
		
		if time_step <= time_lag
			
			for n = 1:length(model_params1)
				model_param1_str = param2str(model_params1(n));
				
				for o = 1:length(model_params2)
					model_param2_str = param2str(model_params2(n));
					
					
					if strcmp(lower(method), 'kraskov') | strcmp(lower(method), 'gaussian');
						
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
								micro_variable_names{k}));
						end
						
						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathout_data_sim_time_series network '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
								macro_variable_names{k}));
						end

					else
						
						bin_number_str = num2str(bin_number);
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							micro_variables(1).([disc_method bin_number '_' micro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' disc_method bin_number '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
								[disc_method bin_number '_' micro_variable_names{k}]));
						end
						
						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							macro_variables(1).([disc_method bin_number '_' macro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' disc_method bin_number '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
								[disc_method bin_number '_' macro_variable_names{k}]));
						end

					end
					
					% calculate dynamical independence (DD)
					fprintf('get all DD - loop indices: time_step: %d, model_param1: %d, model_param2: %d\n', a, n, o);
											
					% get DD for all combinations of micro and top-level macro variables
					DD = get_DD(micro_variables, macro_variables, method, time_lag, time_step, kraskov_param);
					
					%fieldnames_DD = fieldnames(DD);
					
					% big_struct(1,y).results(1,1).([micro_variable_names{u} '_' macro_variable_names{w}]) = [];
					
					for f = 1:length(fieldnames(input_struct(1,1).results))
						output_struct(1,a).results(1,1).([fieldnames_results{f}])(n,o) = {DD.([fieldnames_results{f}])};
					end
					
					clear DD;
					clear macro_variables;
					clear micro_variables;
				end
			end
			
			output_struct(1,a).measure		= measure;
			output_struct(1,a).method		= method;
			output_struct(1,a).time_lag		= time_lag;
			output_struct(1,a).time_length	= time_length;
			output_struct(1,a).kraskov_param	= kraskov_param;
			output_struct(1,a).disc_method	= disc_method;
			output_struct(1,a).bin_number		= bin_number;
			output_struct(1,a).time_step		= time_step;
			
		else
			continue
		end
	end
end 
								
							
