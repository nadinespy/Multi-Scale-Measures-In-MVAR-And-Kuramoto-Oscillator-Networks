function get_all_quant_variables(network, measure_params, model_sim_params, micro_variable_names, ...
		macro_variable_names, pathin_sim_time_series);

	% Function description: get_all_quant_variables() takes as inputs, amongst others, 
	% arrays with measure parameter values common to all emergence calculations 
	% ([measure_params], required), model parameters used for simulation 
	% ([model_sim_params], required).
	% It then loops over time-lengths, model parameters, and number of bins. 
	
	% Inputs:	
	%
	% Required:	network				char array
	%		model_sim_params			1x1 struct with fields 
	%							model_params1, model_params2, 
	%							each of which contain 
	%							float arrays, and other 
	%							possible fields which are not
	%							used here
	%
	%		measure_params			1x1 struct with fields
	%							'measures', 'methods',
	%							'time_lags', 'time_lengths',
	%							'kraskov_params', 'disc_methods',
	%							' bins'
	%
	%							'measures': cell array with chars
	%							'methods': cell array with chars
	%							'time_lags': double array
	%							'time_lengths': double array
	%							'kraskov_params': double array
	%							'disc_methods': cell arrays with chars
	%							'bins': int array
	%
	%		micro_variable_names		cell array with chars
	%		macro_variable_names		cell array with chars
	%		pathin_sim_time_series		1x1 struct with field 
	%							indicating path to input, 
	%							and char array as value
	%
	% Outputs: quantilized time-series of variables specified in 			 
	%	     [micro_variable_names] and [macro_variable_names].							
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'model_sim_params', @isstruct);
	addRequired(p,'measure_params', @isstruct);
	addRequired(p,'micro_variable_names', @iscell);
	addRequired(p,'macro_variable_names', @iscell);
	addRequired(p,'pathin_sim_time_series', @isstruct);
	
	parse(p, network, model_sim_params, measure_params, ...
		micro_variable_names, macro_variable_names, ...
		pathin_sim_time_series);
	
	network				= p.Results.network;
	model_sim_params			= p.Results.model_sim_params;
	measure_params			= p.Results.measure_params;
	micro_variable_names		= p.Results.micro_variable_names;
	macro_variable_names		= p.Results.macro_variable_names;
	pathin_sim_time_series		= p.Results.pathin_sim_time_series;
	
	% extract cell arrays from structs
	fieldnames_pathin			= fieldnames(pathin_sim_time_series);
	pathin				= getfield(pathin_sim_time_series, fieldnames_pathin{1});
	
	model_sim_params_fieldnames	= fieldnames(model_sim_params);
	model_params1			= model_sim_params.(model_sim_params_fieldnames{1});
	model_params2			= model_sim_params.(model_sim_params_fieldnames{2});
	
	time_lengths			= measure_params.time_lengths;
	bins					= measure_params.bins;
	

	
	for q = 1:length(time_lengths);
		time_length_str = num2str(time_lengths(q));
		
		rng(1);
		for i = 1:length(model_params1);
			model_param1_str = param2str(model_params1(i));
			
			for j = 1:length(model_params2)
				model_param2_str = param2str(model_params2(j));
				
				for o = 1:length(bins);
					bin_number = bins(o);
					bin_number_str = num2str(bin_number);
				
					fprintf('get_all_quant_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d, measure_param5: %d\n', q, i, j, o);
					
					% load all micro variables into one struct 'micro_variables'
					for k = 1:length(micro_variable_names)
						micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathin network '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
							micro_variable_names{k}));
					end
						
					% load all macro variables into one struct 'macro_variables'
					for k = 1:length(macro_variable_names)
						macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathin network '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
							macro_variable_names{k}));
					end

					% quantilize all micro variables, and store them in another struct 'quant_micro_variables'
					for k = 1:length(micro_variable_names)
						[quantiles_variable, quantilized_variable] = get_quant_variables(micro_variables(1).(micro_variable_names{k}), bin_number);
						quant_micro_variables(1).(['quant' bin_number_str '_' micro_variable_names{k}]) = quantilized_variable;
					end
					
					% quantilize all macro variables, and store them in another struct 'quant_macro_variables'
					for k = 1:length(macro_variable_names)
						[quantiles_variable, quantilized_variable] = get_quant_variables(macro_variables(1).(macro_variable_names{k}), bin_number);
						quant_macro_variables(1).(['quant' bin_number_str '_' macro_variable_names{k}]) = quantilized_variable;
					end
					
					% get all fieldnames from 'quant_macro_variables' & 'quant_macro_variables'
					quant_micro_variable_names = fieldnames(quant_micro_variables);
					quant_macro_variable_names = fieldnames(quant_macro_variables);
					
					% save all quantilized micro variables from 'quant_micro_variables' as single variables
					% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints + number of quantiles
					for k = 1:length(quant_micro_variable_names);
						% save single element from struct
						save([pathin network '_' quant_micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
							'-struct', 'quant_micro_variables', quant_micro_variable_names{k});
					end
					
					% save all quantilized macro variables from 'quant_macro_variables' as single variables
					% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints + number of quantiles
					for k = 1:length(quant_macro_variable_names);
						% save single element from struct
						save([pathin network '_' quant_macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
							'-struct', 'quant_macro_variables', quant_macro_variable_names{k});
					end
					
				end
			end
		end
	end
	
end
