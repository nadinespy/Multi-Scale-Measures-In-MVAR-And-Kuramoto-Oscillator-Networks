function output_struct = get_all_phiidCE_DC_CD(network, ...
		model_params1, ...
		model_params2, ...
		time_length, ...
		measure, ...
		method, ...
		time_lag, ...
		red_funcs, ...
		micro_variable_names, ...
		input_struct, ...
		pathin_sim_time_series, ...
		varargin)
	
	% Function description: get_all_ShannonCE_DC_CD() takes as inputs scalar parameter 
	% values common to all emergence calculations, and arrays with parameter values
	% specific to PhiID-CE/DC/CD. It then loops overPhiID-CE/DE/CD-specific parameters.
	
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
	%		red_funcs				string array
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
	%							but with PhiID-CE-, or 
	%							PhiID-CD-, or PhiID-DC 
	%							values (as opposed to zeros)
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
	addRequired(p,'red_funcs', @iscell);
	addRequired(p,'micro_variable_names', @iscell);
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
		measure, method, time_lag, red_funcs, micro_variable_names, ...
		input_struct, pathin_sim_time_series, varargin{:});

	network				= p.Results.network;
	model_params1			= p.Results.model_params1;
	model_params2			= p.Results.model_params2;
	measure				= p.Results.measure;
	method				= p.Results.method;
	time_length				= p.Results.time_length;
	time_lag				= p.Results.time_lag;
	red_funcs				= p.Results.red_funcs;
	micro_variable_names		= p.Results.micro_variable_names;
	input_struct			= p.Results.input_struct;
	pathin_sim_time_series		= p.Results.pathin_sim_time_series;
	kraskov_param			= p.Results.kraskov_param;
	disc_method				= p.Results.disc_method;
	bin_number				= p.Results.bin_number;
	
	
	output_struct = input_struct;
	fieldnames_results = fieldnames(input_struct(1,1).results(1,1));
	
	fieldnames_pathin = fieldnames(pathin_sim_time_series);
	pathin = getfield(pathin_sim_time_series, fieldnames_pathin{1});
	
	time_length_str = num2str(time_length);
	bin_number_str = num2str(bin_number);
	
	% loop over time_steps, model_params1, model_params2
		
	for a = 1:length(red_funcs)
		red_func = red_funcs{a};
	
		for n = 1:length(model_params1)
			model_param1_str = param2str(model_params1(n));
		
			for o = 1:length(model_params2)
				model_param2_str = param2str(model_params2(n));
			
				if strcmp(lower(method), 'kraskov') | strcmp(lower(method), 'gaussian');
					
					% load all micro variables into one struct 'micro_variables'
					for k = 1:length(micro_variable_names)
						micro_variables(1).(micro_variable_names{k}) = struct2array(load([pathin network '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
							micro_variable_names{k}));
					end
					
				else
					
					bin_number_str = num2str(bin_number);
					% load all micro variables into one struct 'micro_variables'
					for k = 1:length(micro_variable_names)
						micro_variables(1).([disc_method bin_number_str '_' micro_variable_names{k}]) = struct2array(load([pathin network '_' disc_method bin_number_str '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
							[disc_method bin_number_str '_' micro_variable_names{k}]));
					end
					
				end
			
				% calculate PhiID-CE, PhiID-DC, or PhiID-CD
				fprintf('get all phiid-CE/DC/CD - loop indices: time_step: %d, model_param1: %d, model_param2: %d\n', a, n, o);
			
				% get PhiID-CE, PhiID-DC, or PhiID-CD for all combinations of micro and top-level macro variables: 
				% all_phiidCE_DC_CD is an array with structs where each struct corresponds to one PhiID-measure specified in [measure],
				% in the same order
				all_phiidCE_DC_CD = get_phiidCE_DC_CD(micro_variables, measure, method, time_lag, red_func, ...
					'kraskov_param', kraskov_param);
						
				% big_struct(1,y).results(1,1).([micro_variable_names{u} '_' macro_variable_names{w}]) = [];
				for h = 1:length(all_phiidCE_DC_CD)
					for f = 1:length(fieldnames(input_struct(1,1).results))
						output_struct(1,a).results(1,1).([fieldnames_results{f}])(n,o) = {all_phiidCE_DC_CD{h}.([fieldnames_results{f}])};
						output_struct(1,a).measure						   = measure{h};
					end
				end
				
				clear phiidMeasure;
				clear macro_variables;
				clear micro_variables;
		end
	end
			
	output_struct(1,a).method		= method;
	output_struct(1,a).time_lag		= time_lag;
	output_struct(1,a).time_length	= time_length;
	output_struct(1,a).kraskov_param	= kraskov_param;
	output_struct(1,a).disc_method	= disc_method;
	output_struct(1,a).bin_number		= bin_number;

end 
								
							
