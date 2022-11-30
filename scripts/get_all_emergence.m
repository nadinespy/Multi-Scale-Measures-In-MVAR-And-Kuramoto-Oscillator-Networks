function emergence_struct = get_all_emergence(network, model_calc_params, ...
		measure_params, micro_variable_names, macro_variable_names, struct_prefix, ...
		pathin_sim_time_series, pathout_emergence, varargin)
% get_all_emergence() - calculates one or more measures of emergence for one or more
% models for one or more micro and macro variables, respectively. 
% 
% Takes as inputs arrays with measure parameter values common to all emergence measures 
% ([measure_params], required), model parameters ([model_calc_params], required), and 
% arrays with measure parameter values specific to measures ([varargin]). It loops over 
% measure parameter values common to all emergence measures, calling measure-specific 
% functions and looping over measure-specific parameters which then loop over model 
% parameters.
% 
% Example: emergence_results = get_all_emergence(network, model_calc_params, ...
%		measure_params, micro_variable_names, macro_variable_names, ...
%		struct_prefix, pathin, pathout, 'measure_params_phiid_ce', ...
%		measure_params_phiid_ce, 'measure_params_dd', measure_params_dd);
%
% Inputs - required:	
%    network                -             char array
%    model_calc_params      -             1x1 struct with fields 
%                                         model_params1, model_params2, 
%                                         each of which contain 
%                                         float arrays
%
%    measure_params         -             1x1 struct with fields
%                                         'measures', 'methods',
%                                         'time_lags', 'time_lengths',
%                                         'kraskov_params', 'disc_methods',
%                                         ' bins'
%
%                                         'measures': cell array with chars
%                                         'methods': cell array with chars
%                                         'time_lags': double array
%                                         'time_lengths': double array
%                                         'kraskov_params': double array
%                                         'disc_methods': cell arrays with chars
%                                         'bins': int array
%
%   micro_variable_names    -             cell array with chars
%   macro_variable_names    -             cell array with chars
%   struct_prefix             -             char array
%   pathin_sim_time_series  -             1x1 struct with field 
%                                         indicating path to input, 
%                                         and char array as value
%   pathout_emergence       -             1x1 struct with fields 
%                                         indicating paths to results, e. g., 
%                                         'pathout_data_shannon_ce',
%                                         'pathout_data_phiid_ce', or
%                                         'pathout_data_dd', each of
%                                         which has char array as value
%
% Inputs - optional: 
%   measure_params_dd       -             1x1 struct with fields 
%                                         'time_steps': double array
%   measure_params_phiid_ce -             1x1 struct with fields 
%                                         'red_funcs': cell array with chars
%
% Outputs: 
%   emergence_results       -             1xN struct with fields: 
%                                         'time_length',
%                                         'measure', 
%                                         'method', 
%                                         'time_lag', 
%                                         'disc_method', 
%                                         'bin_number', 
%                                         'kraskov_param',
%                                         'time_step',
%                                         'red_func',
%                                         'results'
%
%                                         where 'results' contains
%                                         fieldnames according to
%                                         micro-macro combinations,
%                                         each of which is valued with 
%                                         a table with emergence results,
%                                         with model_params1
%                                         and model_params2 as rows/
%                                         columns; all other fields
%                                         contain one value from the 
%                                         structs' variables 
%                                         described above; N is the 
%                                         total number of parameter 
%                                         combinations (excluding 
%                                         model parameters) for which 
%                                         emergence is calculated, 
%                                         e. g., N = number of measures * 
%                                                    number of methods * 
%                                                    number of time-lags * 
%                                                    number of time-lengths * 
%                                                    number of K-nearest 
%                                                    neighbours							     
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'model_calc_params', @isstruct);
	addRequired(p,'measure_params', @isstruct);
	addRequired(p,'micro_variable_names', @iscell);
	addRequired(p,'macro_variable_names', @iscell);
	addRequired(p,'struct_prefix', @ischar);
	addRequired(p,'pathin_sim_time_series', @isstruct);
	addRequired(p,'pathout_emergence', @isstruct);
	
	% optional name-value pair variables: 
	default_measure_params_dd = struct('time_steps', [1]);
	addParameter(p,'measure_params_dd', default_measure_params_dd, ...
		@isstruct);
	default_measure_params_phiid_ce = struct('red_funcs', {'MMI', 'CCS'});
	addParameter(p,'measure_params_phiid_ce', ...
		default_measure_params_phiid_ce, @isstruct);
	
	parse(p, network, model_calc_params, measure_params, ...
		micro_variable_names, macro_variable_names, struct_prefix, ...
		pathin_sim_time_series, pathout_emergence, varargin{:});
	
	network				= p.Results.network;
	model_calc_params			= p.Results.model_calc_params;
	measure_params			= p.Results.measure_params;
	micro_variable_names		= p.Results.micro_variable_names;
	macro_variable_names		= p.Results.macro_variable_names;
	struct_prefix			= p.Results.struct_prefix;
	pathin_sim_time_series		= p.Results.pathin_sim_time_series;
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
	
	% put PhiID-measures into separate cell array in [measures] for get_all_phiidCE_DC_CD()
	phiidMeasures = {};
	
	z = 1;
	while (length(measures) ~= 0) 
		
		if strcmp(measures{z}, 'phiidCE') || strcmp(measures{z}, 'phiidDC') || strcmp(measures{z}, 'phiidCD')
			phiidMeasure = measures{z};
			phiidMeasures = [phiidMeasures, measures{z}];
			measures(find(strcmp(measures, phiidMeasure))) = [];
		else 
			z = z+1;
		end 
		
		if z == length(measures)
			break
		end
		
	end 
	measures = [measures, {phiidMeasures}];
	
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
		
	% initialise master structure where to store all results
	emergence_struct = struct('time_length', [], 'measure', [], 'method', [], ...
		'time_lag', [], 'disc_method', [], 'bin_number', [], ...
		'kraskov_param', [], 'time_step', [], 'red_func', [], ...
		'results', struct([]));

	% placeholder structs to fill in emergence results for measure specific parameters; placeholder
	% struct for PhiID-measures must have space not only for red_funcs, but also for different
	% PhiID-measures
	temp_emergence_struct_DD(1,1:length(time_steps))					= emergence_struct;
	temp_emergence_struct_phiid(1,1:length(red_funcs)*length(phiidMeasures))	= emergence_struct;
	temp_emergence_struct										= emergence_struct;
	
	% run the big loop over time lengths, measures, time lags, methods, 
	% discretization methods & number of bins (for method 'discrete'), 
	% K-nearest neighbours (for method 'kraskov'), time steps (in DD), 
	% and reduncancy functions (in PhiID-CE);
	% emergence_struct is appended with emergence results for each 
	% parameter combination (excluding model parameters which are stored 
	% in [emergence_struct.results.('micro_variable_name_macro_variable_name')])
	
	for j = 1:length(time_lengths)
		time_length = time_lengths(j);
		
		for q = 1:length(measures);
			measure = measures{q};
			
			for z = 1:length(methods);
				method = methods{z};
				
				for i = 1:length(time_lags);
					time_lag = time_lags(i);
					
					% ------------------------------------------------------
					% method 'Discrete'
					% ------------------------------------------------------
					if strcmp(lower(method), 'discrete');
						
						for c = 1:length(disc_methods)
							disc_method = disc_methods{c};
							
							for l = 1:length(bins)
								bin_number = bins(l);
								bin_number_str = num2str(bin_number);
								
								fprintf('get all emergence - loop indices: time_length: %d, measure: %d, method: %d, time_lag: %d, disc_method: %d, bin_number: %d\n', j, q, z, i, c, l);
								
								% ------------------------------------------------------
								% PhiID-CE/DC/CD
								% ------------------------------------------------------
								if iscell(measure) == true
									
									if any(strcmp(measure, 'phiidCE')) || any(strcmp(measure, 'phiidDC')) || any(strcmp(measure, 'phiidCD'));
									
										all_phiidCE_DC_CD = get_all_phiidCE_DC_CD(network, ...
											model_params1, ...
											model_params2, ...
											time_length, ...
											measure, ...
											method, ...
											time_lag, ...
											red_funcs, ...
											micro_variable_names, ...
											temp_emergence_struct_phiid(1,1:length(red_funcs)*length(phiidMeasures)), ...
											pathin_sim_time_series, ...
											'disc_method', disc_method, ...
											'bin_number', bin_number);
										
										emergence_struct = [emergence_struct, all_phiidCE_DC_CD];
										
									end 
									
								% ------------------------------------------------------
								% Shannon-CE, Shannon-DC, or Shannon-CD
								% ------------------------------------------------------
								
								elseif strcmp(measure, 'shannonCE') || strcmp(measure, 'shannonDC') || strcmp(measure, 'shannonCD');
									
									all_shannonCE_DC_CD = get_all_shannonCE_DC_CD(network, ...
										model_params1, ...
										model_params2, ...
										time_length, ...
										measure, ...
										method, ...
										time_lag, ...
										micro_variable_names, ...
										macro_variable_names, ...
										temp_emergence_struct, ...
										pathin_sim_time_series, ...
										'disc_method', disc_method, ...
										'bin_number', bin_number);
									
									emergence_struct = [emergence_struct, all_shannonCE_DC_CD];
									
								% ------------------------------------------------------
								% Dynamical Dependence
								% ------------------------------------------------------
								elseif strcmp(measure, 'DD');
									
									% loop over time_steps, model_params1, and model_params2, and append
									% results to emergence_struct
									all_DD = get_all_DD(network, ...
										model_params1, ...
										model_params2, ...
										time_length, ...
										measure, ...
										method, ...
										time_lag, ...
										time_steps, ...
										micro_variable_names, ...
										macro_variable_names, ...
										temp_emergence_struct_DD(1,1:length(time_steps)), ...
										pathin_sim_time_series, ...
										'disc_method', disc_method, ...
										'bin_number', bin_number);
									
									emergence_struct = [emergence_struct, all_DD];
								
								end
							end
						end
						
					% ------------------------------------------------------
					% method 'Kraskov'
					% ------------------------------------------------------
					elseif strcmp(lower(method), 'kraskov');
						
						for h = 1:length(kraskov_params);
							kraskov_param = kraskov_params(h);
							fprintf('get all emergence - loop indices: time_length: %d, measure: %d, method: %d, time_lag: %d, kraskov_param: %d\n', j, q, z, i, h);
							
							% ------------------------------------------------------
							% PhiID-CE/DC/CD
							% ------------------------------------------------------
							if iscell(measure) == true
									
								if any(strcmp(measure, 'phiidCE')) || any(strcmp(measure, 'phiidDC')) || any(strcmp(measure, 'phiidCD'));
									
									all_phiidCE_DC_CD = get_all_phiidCE_DC_CD(network, ...
										model_params1, ...
										model_params2, ...
										time_length, ...
										measure, ...
										method, ...
										time_lag, ...
										red_funcs, ...
										micro_variable_names, ...
										temp_emergence_struct_phiid(1,1:length(red_funcs)*length(phiidMeasures)), ...
										pathin_sim_time_series, ...
										'kraskov_param', kraskov_param);
									
									emergence_struct = [emergence_struct, all_phiidCE_DC_CD];
									
								end
								
							% ------------------------------------------------------
							% Shannon-CE, Shannon-DC, or Shannon-CD
							% ------------------------------------------------------
							elseif strcmp(measure, 'shannonCE') || strcmp(measure, 'shannonDC') || strcmp(measure, 'shannonCD');
								
								all_shannonCE_DC_CD = get_all_shannonCE_DC_CD(network, ...
									model_params1, ...
									model_params2, ...
									time_length, ...
									measure, ...
									method, ...
									time_lag, ...
									micro_variable_names, ...
									macro_variable_names, ...
									temp_emergence_struct, ...
									pathin_sim_time_series, ...
									'kraskov_param', kraskov_param);
								
								emergence_struct = [emergence_struct, all_shannonCE_DC_CD];
								
							% ------------------------------------------------------
							% Dynamical Dependence
							% ------------------------------------------------------
							elseif strcmp(measure, 'DD');
								
								% loop over time_steps, model_params1, and model_params2, and append
								% results to emergence_struct
								all_DD = get_all_DD(network, ...
									model_params1, ...
									model_params2, ...
									time_length, ...
									measure, ...
									method, ...
									time_lag, ...
									time_steps, ...
									micro_variable_names, ...
									macro_variable_names, ...
									temp_emergence_struct_DD(1,1:length(time_steps)), ...
									pathin_sim_time_series, ...
									'kraskov_param', kraskov_param);
								
								emergence_struct = [emergence_struct, all_DD];

							end
						end
						
					% ------------------------------------------------------
					% method 'Gaussian'
					% ------------------------------------------------------
					elseif strcmp(lower(method), 'gaussian')
						
						fprintf('get all emergence - loop indices: time_length: %d, measure: %d, method: %d, time_lag: %d, ', j, q, z, i);
						
						% ------------------------------------------------------
						% PhiID-CE/DC/CD
						% ------------------------------------------------------
						if iscell(measure) == true
							
							if any(strcmp(measure, 'phiidCE')) || any(strcmp(measure, 'phiidDC')) || any(strcmp(measure, 'phiidCD'));
								
								all_phiidCE_DC_CD = get_all_phiidCE_DC_CD(network, ...
									model_params1, ...
									model_params2, ...
									time_length, ...
									measure, ...
									method, ...
									time_lag, ...
									red_funcs, ...
									micro_variable_names, ...
									temp_emergence_struct_phiid(1,1:length(red_funcs)*length(phiidMeasures)), ...
									pathin_sim_time_series);
								
								emergence_struct = [emergence_struct, all_phiidCE_DC_CD];
								
							end
							
						% ------------------------------------------------------
						% Shannon-CE, Shannon-DC, or Shannon-CD
						% ------------------------------------------------------
						elseif strcmp(measure, 'shannonCE') || strcmp(measure, 'shannonDC') || strcmp(measure, 'shannonCD');
							
							all_shannonCE_DC_CD = get_all_shannonCE_DC_CD(network, ...
								model_params1, ...
								model_params2, ...
								time_length, ...
								measure, ...
								method, ...
								time_lag, ...
								micro_variable_names, ...
								macro_variable_names, ...
								temp_emergence_struct, ...
								pathin_sim_time_series);
							
							emergence_struct = [emergence_struct, all_shannonCE_DC_CD];
							
						% ------------------------------------------------------
						% Dynamical Dependence
						% ------------------------------------------------------
						elseif strcmp(measure, 'DD')
							
							% loop over time_steps, model_params1, and model_params2, and append
							% results to emergence_struct
							
							all_DD = get_all_DD(network, ...
								model_params1, ...
								model_params2, ...
								time_length, ...
								measure, ...
								method, ...
								time_lag, ...
								time_steps, ...
								micro_variable_names, ...
								macro_variable_names, ...
								temp_emergence_struct_DD(1,1:length(time_steps)), ...
								pathin_sim_time_series);
							
							emergence_struct = [emergence_struct, all_DD];
						end
					end
				end
			end
		end
	end
	
	emergence_struct = emergence_struct(1,2:end);
end	
