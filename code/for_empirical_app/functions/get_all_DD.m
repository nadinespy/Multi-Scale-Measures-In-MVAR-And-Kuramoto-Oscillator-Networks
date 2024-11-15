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
% get_all_DD() calculates Dynamical Dependence (DD), for all micro-macro 
% variable combinations, all DD-specific measure parameters, and all 
% model parameters.
%
% Takes as inputs the measure ([measure]), method ([method]), scalar measure 
% parameter values common to all emergence calculations ([time_length], [time_lag]),
% model parameter arrays ([model_params1], [model_params2]), DD-specific measure 
% parameters [time_steps], and names of the micro and macro variables 
% ([micro_variable_names], [macro_variable_names]). It loops over the DD-specific
% parameters.
%
% Example 1: output_struct = get_all_DD(network, model_params1, model_params2, ...
%	     time_length, measure, method, time_lag, time_steps, ...
%	     micro_variable_names, macro_variable_names, input_struct, ...
%          pathin_sim_time_series, 'kraskov_param', kraskov_param) 
%
% Example 2: output_struct = get_all_DD(network, model_params1, model_params2, ...
%	     time_length, measure, method, time_lag, time_steps, ...
%	     micro_variable_names, macro_variable_names, input_struct, ...
%          pathin_sim_time_series, 'disc_method', disc_method, ...
%	     'bin_number', bin_number) 
%	
% INPUTS - required: 
%    input_struct		    -			1x(length(time_steps))struct 
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
% INPUTS - optional: 
%    kraskov_param	    -			double
%    disc_method		    -			char array
%    bin_number		    -			double	
%
% OUTPUTS: 
%    output_struct	    -			same struct as input_struct,
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
					
					fprintf('get all DD - loop indices: time_step: %d, model_param1: %d, model_param2: %d\n', a, n, o);
					
					if strcmp(lower(method), 'kraskov') | strcmp(lower(method), 'gaussian');
						
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							
							micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathin network '_' ...
								micro_variable_names{k} '_' ...
								model_param1_str '_' ...
								model_param2_str '_' ...
								time_length_str '.mat'], ...
								micro_variable_names{k}));
							
							% transpose micro variable, if datapoints are stored in rows, and variables in column
							if size(micro_variables(1).(micro_variable_names{k}),1) > size(micro_variables(1).(micro_variable_names{k}),2)
								micro_variables(1).(micro_variable_names{k}) = micro_variables(1).(micro_variable_names{k})';
							end
						end
						
						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							
							macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathin network '_' ...
								macro_variable_names{k} '_' ...
								model_param1_str '_' ...
								model_param2_str '_' ...
								time_length_str '.mat'], ...
								macro_variable_names{k}));
							
							% transpose macro variable, if datapoints are stored in rows, and variable(s) in column
							if size(macro_variables(1).(macro_variable_names{k}),1) > size(macro_variables(1).(macro_variable_names{k}),2)
								macro_variables(1).(macro_variable_names{k}) = macro_variables(1).(macro_variable_names{k})';
							end
						
						end

					elseif strcmp(lower(method), 'discrete');
						
						n_features_micro = 2;
						n_features_macro = 1;
						
						% load all micro variables into one struct 'micro_variables'
						for k = 1:length(micro_variable_names)
							
							micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathin network '_' ...
								micro_variable_names{k} '_' ...
								model_param1_str '_' ...
								model_param2_str '_' ...
								time_length_str '.mat'], ...
								micro_variable_names{k}));
							
							% do dimensionality reduction using RICA to be able to use JIDT for the discrete case
							try
								reconstruction_ica = rica(micro_variables(1).(micro_variable_names{k})', n_features_micro);
								rica_micro = (micro_variables(1).(micro_variable_names{k})' * reconstruction_ica.TransformWeights)';
								micro_variables(1).(micro_variable_names{k}) = rica_micro;
								
								% discretize
								[quantiles_rica_micro, quantilized_rica_micro] = get_quant_variables(rica_micro, 1);
								micro_variables(1).(micro_variable_names{k}) = quantilized_rica_micro;
							catch 
								warning('dim reduction of micro variables for discrete method in DD has not worked');
								micro_variables(1).(micro_variable_names{k}) = NaN;
							end

						end
						
						% load all macro variables into one struct 'macro_variables'
						for k = 1:length(macro_variable_names)
							
							macro_variables(1).(macro_variable_names{k}) = struct2array(load([pathin network '_' ...
								macro_variable_names{k} '_' ...
								model_param1_str '_' ...
								model_param2_str '_' ...
								time_length_str '.mat'], ...
								macro_variable_names{k}));
							
							% transpose macro variable, if datapoints are stored in rows, and variable(s) in column
							if size(macro_variables(1).(macro_variable_names{k}),1) > size(macro_variables(1).(macro_variable_names{k}),2)
								macro_variables(1).(macro_variable_names{k}) = macro_variables(1).(macro_variable_names{k})';
							end
							
							% if macro variable has more than three dimensions, reduce them to two to be able to use JIDT for discrete case
							if size(macro_variables(1).(macro_variable_names{k}),1) > 1;
								
								try
									reconstruction_ica = rica(macro_variables(1).(macro_variable_names{k})', n_features_macro);
									rica_macro = (macro_variables(1).(macro_variable_names{k})' * reconstruction_ica.TransformWeights)';
									macro_variables(1).(macro_variable_names{k}) = rica_macro;
									
									% discretize
									[quantiles_rica_macro, quantilized_rica_macro] = get_quant_variables(rica_macro, 1);
									macro_variables(1).(macro_variable_names{k}) = quantilized_rica_macro;
								catch
									warning('dim reduction of macro variables for discrete method in DD has not worked');
									macro_variables(1).(macro_variable_names{k}) = NaN;
								end
								
							else 
								% discretize
								[quantiles_macro, quantilized_macro] = get_quant_variables(macro_variables(1).(macro_variable_names{k}), 1);
								macro_variables(1).(macro_variable_names{k}) = quantilized_macro;
							end
							
						end

					end
					
					% calculate dynamical independence (DD)
											
					% get DD for all combinations of micro and top-level macro variables
					DD = get_DD(micro_variables, macro_variables, method, time_lag, time_step, 'kraskov_param', kraskov_param);
					
					fieldnames_DD = fieldnames(DD);
					
					for f = 1:length(fieldnames_results)
						output_struct(1,a).results(1,1).([fieldnames_results{f}])(n,o) = {DD.([fieldnames_DD{f}])};
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
			output_struct(1,a).model		= network;
			
		else
			continue
		end
	end
end 
								
							
