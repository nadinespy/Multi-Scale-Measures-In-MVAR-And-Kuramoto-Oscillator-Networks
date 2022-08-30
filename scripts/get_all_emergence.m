function emergence_struct = get_all_emergence(network, model_calc_params, ...
		measure_params, micro_variable_names, macro_variable_names, variable_name, ...
		pathout_data_sim_time_series, pathout_emergence, varargin)

	% Inputs:	
	%
	% Required:	network				string
	%		model_calc_params			1x1 struct with fields 
	%							model_params1, model_params2, 
	%							each of which contain 
	%							float arrays
	%
	%		measure_params			1x1 struct with fields
	%							'measures', 'methods',
	%							'time_lags', 'time_lengths',
	%							'kraskov_params', 'disc_methods',
	%							' bins'
	%
	%							'measures': cell array with chars
	%							'methods':cell array with chars
	%							'time_lags': int array
	%							'time_lengths': int array
	%							'kraskov_params': int array
	%							'disc_methods': cell arrays with chars
	%							'bins': int array
	%
	%		micro_variable_names		cell array with chars
	%		macro_variable_names		cell array with chars
	%		variable_name			char
	%		pathout_data_sim_time_series	char
	%		pathout_emergence			1x1 struct with fields 
	%							'pathout_data_shannon_ce',
	%							'pathout_data_phiid_ce',
	%							'pathout_data_dd', each of
	%							which contain a char array
	%
	% Optional: measure_params_dd			1x1 struct with fields 
	%							'time_steps': int array
	%		measure_params_phiid_ce		1x1 struct with fields 
	%							'red_funcs': cell array with chars
	%
	% Outputs: emergence_results			struct with fields 
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
	%							table with emergence results,
	%							with model_params1
	%							and model_params2 as rows/
	%							columns; all oother fields
	%							contain one value from the 
	%							structs' variables 
	%							described above
	
	p = inputParser;
	addRequired(p,'network', @ischar);
	addRequired(p,'model_calc_params', @isstruct);
	addRequired(p,'measure_params', @isstruct);
	addRequired(p,'micro_variable_names', @iscell);
	addRequired(p,'macro_variable_names', @iscell);
	addRequired(p,'variable_name', @ischar);
	addRequired(p,'pathout_data_sim_time_series', @ischar);
	addRequired(p,'pathout_emergence', @isstruct);
	
	default_measure_params_dd = struct('time_steps', [1]);
	addOptional(p,'measure_params_dd', default_measure_params_dd, ...
		@isstruct);
	default_measure_params_phiid_ce = struct('red_funcs', {'MMI', 'CCS'});
	addOptional(p,'measure_params_phiid_ce', ...
		default_measure_params_phiid_ce, @isstruct);
	
	parse(p, network, model_calc_params, measure_params, ...
		micro_variable_names, macro_variable_names, variable_name, ...
		pathout_data_sim_time_series, pathout_emergence, varargin{:});
	
	network				= p.Results.network;
	model_calc_params			= p.Results.model_calc_params;
	measure_params			= p.Results.measure_params;
	micro_variable_names		= p.Results.micro_variable_names;
	macro_variable_names		= p.Results.macro_variable_names;
	variable_name			= p.Results.variable_name;
	pathout_data_sim_time_series	= p.Results.pathout_data_sim_time_series;
	pathout_emergence			= p.Results.pathout_emergence;
	measure_params_dd			= p.Results.measure_params_dd;
	measure_params_phiid_ce		= p.Results.measure_params_phiid_ce;
	
	model_calc_params_fieldnames	= fieldnames(model_calc_params);
	measure_params_fieldnames	= fieldnames(measure_params);
	
	% extract cell arrays from structs
	model_params1	= model_calc_params.(model_calc_params_fieldnames{1});
	model_params2	= model_calc_params.(model_calc_params_fieldnames{2});
	
	measures		= measure_params.measures;
	methods		= measure_params.methods;
	time_lags		= measure_params.time_lags;	
	time_lengths	= measure_params.time_lengths;
	time_steps		= measure_params_dd.time_steps;
	red_funcs		= measure_params_phiid_ce.red_funcs;
	
	if length(measure_params_fieldnames) > 4
		for g = 1:length(measure_params_fieldnames)
			if strcmp(measure_params_fieldnames{g}, 'kraskov_params') == true;
				kraskov_params = measure_params.kraskov_params;
			end
			if strcmp(measure_params_fieldnames{g}, 'disc_methods') == true;
				disc_methods = measure_params.disc_methods;
				bins = measure_params.bins;
			end
		end
	end
	
% 	if exist('measure_params_dd')
% 		time_steps = measure_params_dd.time_steps;
% 	end
% 	
% 	if exist('measure_params_phiid_ce')
% 		red_funcs = measure_params_phiid_ce.red_funcs;
% 	end 

	model_params1_str = {};
	for t = 1:length(model_params1)
		model_params1_str{t} = num2str(model_params1(t));
	end 
		
	for e = 1:length(model_params2)
		model_params2_str{e} = num2str(model_params1(e));
	end 
		
	
	% initialise master structure
	emergence_struct = struct('time_length', [], 'measure', [], 'method', [], ...
		'time_lag', [], 'disc_method', [], 'bin_number', [], ...
		'kraskov_param', [], 'time_step', [], 'red_func', [], ...
		'results', struct([]));
	
	% initialise fields with micro-macro combinations in  
	% [emergence_struct.results] as well as fields for model parameters in 
	% [emergence_struct.results.([micro_variable_names{u} '_' macro_variable_names{w}])]
	
	for u = 1:length(micro_variable_names);
		
		for w = 1:length(macro_variable_names);
			
			emergence_struct(1,1).results(1,1).([micro_variable_names{u} ...
				'_' macro_variable_names{w}]) = ...
				zeros(length(model_params1),length(model_params2));
			
			emergence_struct(1,1).results(1,1).([micro_variable_names{u} ...
				'_' macro_variable_names{w}]) = ...
				array2table(emergence_struct(1,1).results(1,1).([micro_variable_names{u} ...
				'_' macro_variable_names{w}]), 'RowNames', model_params1_str, ...
				'VariableNames', model_params2_str);
			
		end
	end

	% placeholder structs to fill in emergence results for measure specific parameters
	temp_emergence_struct_time_steps(1,1:length(time_steps)) = emergence_struct;
	temp_emergence_struct_red_funcs(1,1:length(red_funcs)) = emergence_struct;
	
	fieldnames_results = fieldnames(emergence_struct(1,1).results(1,1));
	
	% run the big loop over time lengths, measures, time lags, methods, 
	% discretization methods & number of bins (for method 'discrete'), 
	% K-nearest neighbours (for method 'kraskov'), time steps (in DD), 
	% and reduncancy functions (in PhiIDCE)!
	
	for j = 1:length(time_lengths)
		time_length = time_lengths(j);
		
		for q = 1:length(measures);
			measure = measures{q};
			
			for z = 1:length(methods);
				method = methods{z};
				
				for i = 1:length(time_lags);
					time_lag = time_lags(i);
					
					% ------------------------------------------------------
					% method 'discrete'
					% ------------------------------------------------------
					if strcmp(lower(method), 'discrete');
						
						for c = 1:length(disc_methods)
							disc_method = disc_methods{c};
							
							for l = 1:length(bins)
								bin_number = bins(l);
								bin_number_str = num2str(bin_number);
								fprintf('get all emergence - loop indices: time_length: %d, measure: %d, method: %d, time_lag: %d, disc_method: %d, bin_number: %d\n', j, q, z, i, c, l);
								
								% ------------------------------------------------------
								% PhiID CE
								% ------------------------------------------------------
								if strcmp(measure, 'PhiIDCE');
									
									% do calculation for PhiID CE
									
									% ------------------------------------------------------
									% Shannon CE
									% ------------------------------------------------------
								elseif strcmp(measure, 'ShannonCE');
									
									% do calculation for Shannon Ce
									
									% ------------------------------------------------------
									% dynamical dependence
									% ------------------------------------------------------
								elseif strcmp(measure, 'DD');
									
									% include this at some point as optional argument
									kraskov_param = [];
									
									% loop over time_steps, model_params1, and model_params2
									emergence_struct(1,end:length(time_steps)) = get_all_DD(temp_emergence_struct_time_steps(1,1:length(time_steps)), ...
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
									
									emergence_struct(1,end+length(time_steps))
								end
							end
						end
						
						% ------------------------------------------------------
						% method 'kraskov'
						% ------------------------------------------------------
					elseif strcmp(lower(method), 'kraskov');
						
						for h = 1:length(kraskov_params);
							kraskov_param = kraskov_params(h);
							fprintf('get all emergence - loop indices: time_length: %d, measure: %d, method: %d, time_lag: %d, kraskov_param: %d\n', j, q, z, i, h);
							
							% ------------------------------------------------------
							% PhiID CE
							% ------------------------------------------------------
							if strcmp(measure, 'PhiIDCE');
								
								% do calculation for PhiID CE
								
								% ------------------------------------------------------
								% Shannon CE
								% ------------------------------------------------------
							elseif strcmp(measure, 'ShannonCE');
								
								% do calculation for Shannon Ce
								
								% ------------------------------------------------------
								% dynamical dependence
								% ------------------------------------------------------
							elseif strcmp(measure, 'DD');
								
								% include those at some point as optional arguments
								disc_method = [];
								bin_number = [];
								
								% loop over time_steps, model_params1, and model_params2
								emergence_struct(1,end:length(time_steps)) = get_all_DD(temp_emergence_struct_time_steps(1,1:length(time_steps)), ...
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
								
								
							end
							% alternative
							% emergence = measure_function(inputs)
						end
						
						% ------------------------------------------------------
						% method 'gaussian'
						% ------------------------------------------------------
					elseif strcmp(lower(method), 'gaussian')
						
						fprintf('get all emergence - loop indices: time_length: %d, measure: %d, method: %d, time_lag: %d, ', j, q, z, i);
						
						% ------------------------------------------------------
						% PhiID CE
						% ------------------------------------------------------
						if strcmp(measure, 'PhiIDCE')
							
							% do calculation for PhiID CE
							
							% ------------------------------------------------------
							% Shannon CE
							% ------------------------------------------------------
						elseif strcmp(measure, 'ShannonCE')
							
							% do calculation for Shannon Ce
							
							% ------------------------------------------------------
							% dynamical dependence
							% ------------------------------------------------------
						elseif strcmp(measure, 'DD')
							
							% do calculation for DD
						end
						
						% alternative
						% emergence = measure_function(inputs)
						
					end
				end
			end
		end
	end

end	