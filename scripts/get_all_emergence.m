% model name for saving files (type and size of model)
network =				'12km';

%% model parameter specification 

% kuramoto parameter specification
intra_comm_size =			4;					% intra-community size
n_communities =			3;					% number of communities		
A =					linspace(0.08, 0.8, 10);	% vector with different values for A
beta =				linspace(0.08, 0.8, 10);	% vector with different values for noise correlation: use beta values only 
										% up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
										% values of beta; in these cases, sigma chi will be a non-varying 
										% zero macro variable, yielding erroneous values for emergence 
 
%% measure parameter specification
% -------------------------------------------------------------------------
% necessary input arguments

measures =				{'ShannonCE', 'PhiIDCE', 'DD'};			% emergence measures
methods =				{'Gaussian', 'Kraskov', 'Discrete'};		% to be expanded with 'Kraskov' for practCE
time_lags =				[1, 3, 10];							% time-lags
time_lengths =			[2000, 10000];

% -------------------------------------------------------------------------
% optional input arguments (depending on method)

kraskov_params =			[2, 3, 4];							% not yet implemented in loops
disc_methods =			{'quant'}; %, 'even', 'bin'};				% choose discretization method: 'quant' for using quantiles, 
													% 'even' for discretizing into evenly spaced sections of the state space, 
													% 'bin' for binarizing (scripts for latter two not yet modified)
										
bins =				[1, 3, 7];							% number of bins to do discretization for method 'quant' and 'disc'

% -------------------------------------------------------------------------
% input arguments specific to measures

% PhiID-CE
red_funcs =				{'MMI', 'CCS'};
% DD
time_steps =			[1, 3, 10]; 

%% put all parameters into cell structures

% -------------------------------------------------------------------------
% model parameters

% model parameters for simulating kuramoto oscillators; must be in that order
model_sim_params.A			= A ;		 
model_sim_params.beta			= beta;				
model_sim_params.intra_comm_size	= intra_comm_size;
model_sim_params.n_communities	= n_communities;

% model parameters to calculate emergence for; must be in that order
model_calc_params.A			= A ;		 
model_calc_params.beta			= beta;								
														

% -------------------------------------------------------------------------														
% put all common measure parameters common to Shannon CE, PhiID-CE & DD 
% into one cell structure

measure_params.measures = measures;
measure_params.methods = methods;
measure_params.time_lags = time_lags;
measure_params.time_lengths = time_lengths;
measure_params.kraskov_params = kraskov_params;
measure_params.disc_methods = disc_methods;
measure_params.bins = bins;

% -------------------------------------------------------------------------
% put measure parameters specific to DD into one cell structure
measure_params_dd.time_steps = time_steps;

% put measure parameters specific to PhiID-CE into one cell structure
measure_params_phiid_ce.red_funcs = red_funcs;

% -------------------------------------------------------------------------
% pathouts for output
pathout_emergence.pathout_data_shannon_ce = pathout_data_shannon_ce;
pathout_emergence.pathout_data_phiid_ce = pathout_data_phiid_ce;
pathout_emergence.pathout_data_dd = pathout_data_dd;

% pathin
pathout_input_sim_time_series.pathout_data_sim_time_series = pathout_data_sim_time_series;

% -------------------------------------------------------------------------
% group names of variables generated in get_all_variables() into 
% micro and macro variabels
micro_variable_names = {'raw', 'phase', 'sync', 'rica6_phase', ...
	'rica12_phase'};
macro_variable_names = {'mp_sync', 'chi', 'sum_phase', ...
	'sum_rica6_phase', 'sum_rica12_phase'};

% file prefixes to distinguish different struct files, and not overwrite them
variable_name = 'standard';

%% calculate emergence

function emergence_results = get_all_emergence(network, model_params, measure_params, ...
		measure_params_dd, measure_params_phiid_ce, ...
		micro_variable_names, macro_variable_names, variable_name, ...
		pathout_data_sim_time_series, pathout_emergence)

	p = inputParser;
	addRequired(p,'network', @isstring);
	addRequired(p,'model_params', @isstruct);
	addRequired(p,'measure_params', @isstruct);
	addRequired(p,'micro_variable_names', @isstruct);
	addRequired(p,'macro_variable_names', @isstruct);
	addRequired(p,'variable_name', @isstring);
	addRequired(p,'pathout_data_sim_time_series', @isstruct);
	addRequired(p,'pathout_emergence', @isstruct);
	
	default_measure_params_dd = struct();
	addOptional(p,'measure_params_dd', default_measure_params_dd, ...
		@isstruct);
	default_measure_params_phiid_ce = struct();
	addOptional(p,'measure_params_phiid_ce', ...
		default_measure_params_phiid_ce, @isstruct);
	
	parse(p, network, model_params, measure_params, ...
		micro_variable_names, macro_variable_names, variable_name, ...
		pathout_data_sim_time_series, pathout_emergence, varargin{:});
	
	network  =					p.Results.network;
	model_params =				p.Results.model_params;
	measure_params  =				p.Results.measure_params;
	micro_variable_names =			p.Results.micro_variable_names;
	macro_variable_names =			p.Results.macro_variable_names;
	variable_name =				p.Results.variable_name;
	pathout_data_sim_time_series =	p.Results.pathout_data_sim_time_series;
	pathout_emergence =			p.Results.pathout_emergence;
	% measure_params_dd =			p.Results.measure_params_dd;
	% measure_params_phiid_ce =		p.Results.measure_params_phiid_ce;
	
	model_calc_params_fieldnames	= fieldnames(model_calc_params);
	measure_params_fieldnames	= fieldnames(measure_params);
	% measures		= measure_params.(strcmp(measure_params, 'measures'));
	
	model_param1	= model_calc_params.(model_calc_params_fieldnames{1});
	model_param2	= model_calc_params.(model_calc_params_fieldnames{2});
	
	measures		= measure_params.measures;
	methods		= measure_params.methods;
	time_lags		= measure_params.time_lags;	
	time_lengths	= measure_params.time_lengths;
	
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
	
	if exist measure_params_dd
		time_steps = measure_params_dd.time_steps;
	end
	
	if exist measure_params_phiid_ce
		red_funcs = measure_params_phiid_ce.red_funcs;
	end 

	for n = 1:length(model_param1)
		model_param1_str = param2str(model_param1(n));
		
		for o = 1:length(model_param2)
			model_param2_str = param2str(model_param2(n));
			
			for j = 1:length(time_lengths)
				time_length = time_lengths(j);
				time_length_str = num2str(time_length);
				
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
				
				for q = 1:length(measures);
					measure = measures{q};
					% measure_function = @DD;		% define measure function to get
					% rid of redundant loops for measures
					% below
					for z = 1:length(methods);
						method = methods{z};
						
						for i = 1:length(time_lags);
							time_lag = time_lags(i);
							
							if strcmp(lower(method), 'discrete');
								
								for k = 1:length(disc_methods)
									disc_method = disc_methods{k};
									
									for l = 1:length(bins)
										bin = bins(l);
										
										% load all micro variables into one struct 'micro_variables'
										for k = 1:length(micro_variable_names)
											micro_variables(1).([disc_method '_' micro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' disc_method bin '_' micro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
												[disc_method bin '_' micro_variable_names{k}]));
										end
										
										% load all macro variables into one struct 'macro_variables'
										for k = 1:length(macro_variable_names)
											macro_variables(1).([disc_method '_' macro_variable_names{k}]) = struct2array(load([pathout_data_sim_time_series network '_' disc_method bin '_' macro_variable_names{k} '_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
												[disc_method bin '_' macro_variable_names{k}]));
										end
										
										if strcmp(measure, 'PhiIDCE')
											
											% do calculation for PhiID CE
										elseif strcmp(measure, 'ShannonCE')
											
											% do calculation for Shannon Ce
											
											% ---------------------------------------------------------------------------------------------------------------------------------------
											% dynamical dependence
											% ---------------------------------------------------------------------------------------------------------------------------------------
										elseif strcmp(measure, 'DD')
											
											for a = 1:length(time_steps)
												time_step = time_steps(a)
												
												if time_step <= time_lag
													
													% calculate dynamical independence (DD)
%% CONTINUE HERE
													fprintf('get_all_DD - loop indices: time_series_length: %d, measure_param1: %d, model_param1: %d, model_param2: %d, measure_param2: %d, measure_param3: %d, measure_param1_dd: %d\n', q, z, i, j, b, c, d);
													
													% get DD for all combinations of micro and top-level macro variables
													DD = get_DD(micro_variables, macro_variables, method, time_lag, time_step);
													
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
													
													clear DD;
													clear macro_variables;
													clear micro_variables;
												else
													continue
												end
												% do calculation for DD
											end
										end
									end
								end
								
								% alternative
								% emergence = measure_function(inputs)
								
							elseif strcmp(lower(method), 'kraskov')
								
								for k = 1:length(kraskov_params)
									
									for n = 1:length(model_param1)
										
										for o = 1:length(model_param2)
											
											if strcmp(measure, 'PhiIDCE')
												
												% do calculation for PhiID CE
											elseif strcmp(measure, 'ShannonCE')
												
												% do calculation for Shannon Ce
											elseif strcmp(measure, 'DD')
												
												% do calculation for DD
											end
											
											% alternative
											% emergence = measure_function(inputs)
										end
									end
								end
								
							elseif strcmp(lower(method), 'gaussian')
								
								for n = 1:length(model_param1)
									
									for o = 1:length(model_param2)
										
										if strcmp(measure, 'PhiIDCE')
											
											% do calculation for PhiID CE
										elseif strcmp(measure, 'ShannonCE')
											
											% do calculation for Shannon Ce
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
			end
		end
	end
end

								


			
%% saving practical measures for different micro & macro variables

% saved filenames consist of
% network name + type of causal emergence + number of datapoints + time-lag
save([pathout_data_dd network '_all_DD_' dd_variable_name '_' num2str(time_series_length(q)) '_' num2str(measure_param1(z)) '.mat'], ...
	'all_DD');
				
clear all_DD;
	
	
	

