function output_struct = get_all_DD(network, ...
		model_params1, ...
		model_params2, ...
		time_length, ...
		measure, ...
		method, ...
		time_lag, ...
		time_steps, ...
		micro_variable_names, ...
		macro_variable_names, ...
		input_struct, ...
		pathin_sim_time_series, ...
		varargin)
	
	% Function description: get_all_DD() takes as inputs scalar parameter 
	% values common to all emergence calculations, and arrays with 
	% parameter values specific to DD. It then loops over the DD-specific
	% parameters.
	
	% Inputs:	
	%
	% Required: input_struct			1x(length(time_steps))struct 
	%							with fields 
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
	%							each of which is valued with  
	%							a table with model_params1
	%							and model_params2 as rows/
	%							columns, initialized with 
	%							zeros; all other fields
	%							are empty
	% 
	%		model_params1			float array
	%		model_params2			float array
	%		measure				char array
	%		time_length				double
	%		method				char array
	%		time_lag				double
	%		time_steps				double array
	%		micro_variable_names		cell array with chars
	%		macro_variable_names		cell array with chars
	%		network				char array
	%		pathin_sim_time_series		1x1 struct with field 
	%							indicating path to input, 
	%							and char array as value
	%
	% Optional: kraskov_param			double
	%		disc_method				char array
	%		bin_number				double	
	%
	% Outputs: output_struct			same struct as input_struct,
	%							but with DD values (as 
	%							opposed to zeros) in where
	%							where 'results' contains
	%							fieldnames according to
	%							micro-macro combinations,
	%							each of which contain a 
	%							table with emergence results,
	%							with model_params1
	%							and model_params2 as rows/
	%							columns; all other fields
	%							contain one value from the 
	%							structs' variables 
	%							described above;
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'model_params1', @isvector);
	addRequired(p,'model_params2', @isvector);
	addRequired(p,'time_length', @isdouble);
	addRequired(p,'measure', @ischar);
	addRequired(p,'method', @ischar);
	addRequired(p,'time_lag', @isdouble);
	addRequired(p,'time_steps', @isvector);
	addRequired(p,'micro_variable_names', @iscell);
	addRequired(p,'macro_variable_names', @iscell);
	addRequired(p,'input_struct', @isstruct);
	addRequired(p,'pathin_sim_time_series', @isstruct);
	
	% optional name-value pair variables: 
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, @isdouble);
	default_disc_method = 'quant';
	addParameter(p,'disc_method', default_disc_method, @ischar);
	default_bin_number = 3;
	addParameter(p,'bin_number', default_bin_number, @isdouble);
	
	parse(p, network, model_params1, model_params2, time_length, ...
		measure, method, time_lag, time_steps, ...
		micro_variable_names, macro_variable_names, input_struct, ...
		pathin_sim_time_series, varargin{:});

	network				= p.Results.network;
	model_params1			= p.Results.model_params1;
	model_params2			= p.Results.model_params2;
	measure				= p.Results.measure;
	method				= p.Results.method;
	time_length				= p.Results.time_length;
	time_lag				= p.Results.time_lag;
	time_steps				= p.Results.time_steps;
	micro_variable_names		= p.Results.micro_variable_names;
	macro_variable_names		= p.Results.macro_variable_names;
	input_struct			= p.Results.input_struct;
	pathin_sim_time_series		= p.Results.pathin_sim_time_series;
	kraskov_param			= p.Results.kraskov_param;
	disc_method				= p.Results.disc_method;
	bin_number				= p.Results.bin_number;
	
	output_struct = input_struct;
	
	% get row and column names for table in 
	% [emergence_struct.results.('micro_variable_name_macro_variable_name')]
	model_params1_str = {};
	for t = 1:length(model_params1)
		model_params1_str{t} = num2str(model_params1(t));
	end 
		
	for e = 1:length(model_params2)
		model_params2_str{e} = num2str(model_params1(e));
	end 
	
		% initialise fields with 'micro_variable_name_macro_variable_name' in
	% [output_struct.results] as well as tables with model parameters
	% in rows and columns for each
	% [output_struct.results.('micro_variable_name_macro_variable_name')];
	% same for fields with 'micro_variable_name' for macro-independent PhiID
	
	for u = 1:length(micro_variable_names);
		
		for w = 1:length(macro_variable_names);
			
			for z = 1:length(output_struct)
			
				output_struct(1,z).results(1,1).([micro_variable_names{u} ...
					'_' macro_variable_names{w}]) = ...
					zeros(length(model_params1),length(model_params2));
				
				output_struct(1,z).results(1,1).([micro_variable_names{u} ...
					'_' macro_variable_names{w}]) = ...
					array2table(output_struct(1,z).results(1,1).([micro_variable_names{u} ...
					'_' macro_variable_names{w}]), 'RowNames', model_params1_str, ...
					'VariableNames', model_params2_str);
			end
		end
	end
	
	% output_struct = input_struct;
	fieldnames_results = fieldnames(output_struct(1,1).results(1,1));
	
	fieldnames_pathin = fieldnames(pathin_sim_time_series);
	pathin = getfield(pathin_sim_time_series, fieldnames_pathin{1});
	
	time_length_str = num2str(time_length);
	bin_number_str = num2str(bin_number);
	
	% loop over time_steps, model_params1, model_params2
	for a = 1:length(time_steps)
		time_step = time_steps(a);
		
		if time_step <= time_lag
			
			for n = 1:length(model_params1)
				model_param1_str = param2str(model_params1(n));
				
				for o = 1:length(model_params2)
					model_param2_str = param2str(model_params2(o));
					
					
					if strcmp(lower(method), 'kraskov') | strcmp(lower(method), 'gaussian');
						
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

					elseif strcmp(lower(method), 'discrete');
						
						n_features = 2;
						
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathin network '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
								micro_variable_names{k}));
							
							% do dimensionality reduction using RICA, then binarize to be able to use JIDT for the discrete case
							reconstruction_ica = rica(micro_variables(1).(micro_variable_names{k})', n_features);
							rica_micro = (micro_variables(1).(micro_variable_names{k})' * reconstruction_ica.TransformWeights)';
							[quantiles_rica_micro, quantilized_rica_micro] = get_quant_variables(rica_micro, 1);
							micro_variables(1).(micro_variable_names{k}) = quantilized_rica_micro;
						end
						
						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathin network '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
								macro_variable_names{k}));
							[quantiles_macro, quantilized_macro] = get_quant_variables(macro_variables(1).(macro_variable_names{k}), 1);
							macro_variables(1).(macro_variable_names{k}) = quantilized_macro;
						end

					end
					
					% calculate dynamical independence (DD)
					fprintf('get all DD - loop indices: time_step: %d, model_param1: %d, model_param2: %d\n', a, n, o);
											
					% get DD for all combinations of micro and top-level macro variables
					DD = get_DD(micro_variables, macro_variables, method, time_lag, time_step, 'kraskov_param', kraskov_param);
					
					%fieldnames_DD = fieldnames(DD);
					
					% big_struct(1,y).results(1,1).([micro_variable_names{u} '_' macro_variable_names{w}]) = [];
					
					for f = 1:length(fieldnames_results)
						output_struct(1,a).results(1,1).([fieldnames_results{f}])(n,o) = {DD.([fieldnames_results{f}])};
					end
					
					clear DD;
					clear macro_variables;
					clear micro_variables;
				end
			end
			
			if strcmp(lower(method), 'discrete')
				
				output_struct(1,a).disc_method	= disc_method;
				output_struct(1,a).bin_number		= bin_number;
				
			elseif strcmp(lower(method), 'kraskov')
				
				output_struct(1,a).kraskov_param	= kraskov_param;
			end 
				
			output_struct(1,a).measure		= measure;
			output_struct(1,a).method		= method;
			output_struct(1,a).time_lag		= time_lag;
			output_struct(1,a).time_length	= time_length;
			output_struct(1,a).time_step		= time_step;
			
		else
			continue
		end
	end
end 
								
							
