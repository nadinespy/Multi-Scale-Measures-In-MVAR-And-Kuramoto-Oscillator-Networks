function all_shannonCE_DC_CD = get_shannonCE_DC_CD(micro_variables, macro_variables, ...
		measure, method, time_lag, varargin);
% get_ShannonCE_DC_CD() calculates Shannon-CE, or Shannon-DC, or 
% Shannon-CD for any micro and macro variable specified 
% in [micro_variables] and [macro_variables], and the parameters given; 
% 
% Example: all_shannonCE_DC_CD = get_shannonCE_DC_CD(micro_variables, ...
%		macro_variables, measure, method, time_lag, 'kraskow_param', ...
%		kraskov_param)
%
% INPUTS - required:	
%    micro_variables	    -			1x1 struct with micro 
%							variable names in fields, 
%							and micro variable 
%							time-series as values
%    macro_variables	    -			1x1 struct with macro 
%							variable names in fields, 
%							and macro variable 
%							time-series as values
%    measure		    -			character array
%    method			    -			character array
%    time_lag		    -			double
%
% INPUTS - optional: 
%    kraskov_param	    -			double
%
% OUTPUTS:	
%    all_ShannonCE_DC_CD    -			1x1 struct with with micro- 
%							macro combinations in  
%							fields, and and Shannon-CE, or 
%							Shannon-DC, or Shannon-CD as 
%							values

	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'micro_variables', @isstruct);
	addRequired(p,'macro_variables', @isstruct);
	addRequired(p,'measure', @ischar);
	addRequired(p,'method', @ischar);
	addRequired(p,'time_lag', @isdouble);
	
	% optional name-value pair variables: 
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, @isdouble);
	
	parse(p, micro_variables, macro_variables, measure, method, ...
		time_lag, varargin{:});
	
	measure		= p.Results.measure;
	method		= p.Results.method;
	time_lag		= p.Results.time_lag;
	micro_variables	= p.Results.micro_variables;
	macro_variables	= p.Results.macro_variables;
	kraskov_param	= p.Results.kraskov_param;

	% get number of micro and macro variables, respectively
	n_micro_variables = length(fieldnames(micro_variables)); 
	n_macro_variables = length(fieldnames(macro_variables)); 
	
	% get names of micro and macro variables, respectively
	fieldnames_macro = fieldnames(macro_variables);
	fieldnames_micro = fieldnames(micro_variables);
							
	% loop over all micro and macro variables and calculate practical CE
	for i = 1:n_micro_variables;
		micro = micro_variables.(fieldnames_micro{i});
		
		for j = 1:n_macro_variables;
			macro = macro_variables.(fieldnames_macro{j});
			
			if strcmp(lower(method), 'kraskov')
				
				if strcmp(measure, 'shannonCE');
					try
						shannonMeasure = EmergencePsi(micro', macro', time_lag, method, 'kraskov_param', kraskov_param);
					catch 
						shannonMeasure = NaN;
					end
					
				elseif strcmp(measure, 'shannonDC');
					try
						shannonMeasure = EmergenceDelta(micro', macro', time_lag, method, 'kraskov_param', kraskov_param);
					catch
						shannonMeasure = NaN;
					end
						
				elseif strcmp(measure, 'shannonCD');
					try
						shannonMeasure = EmergenceGamma(micro', macro', time_lag, method, 'kraskov_param', kraskov_param);
					catch
						shannonMeasure = NaN;
					end
				end
				
			else
				
				if strcmp(measure, 'shannonCE')
					try
						shannonMeasure = EmergencePsi(micro', macro', time_lag, method);
					catch
						shannonMeasure = NaN;
					end
					
				elseif strcmp(measure, 'shannonDC')
					try
						shannonMeasure = EmergenceDelta(micro', macro', time_lag, method);
					catch
						shannonMeasure = NaN;
					end
					
				elseif strcmp(measure, 'shannonCD')
					try
						shannonMeasure = EmergenceGamma(micro', macro', time_lag, method);
					catch
						shannonMeasure = NaN;
					end
				end
				
			end
				
			% calculate practical CE, DC, and CE for a given micro-macro combination, and store 
			% results in struct
				
			all_shannonCE_DC_CD.([fieldnames_micro{i} '_' fieldnames_macro{j}]) = shannonMeasure;
		end 
	end 
end 