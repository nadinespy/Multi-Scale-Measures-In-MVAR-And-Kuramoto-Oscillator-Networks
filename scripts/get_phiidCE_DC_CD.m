function all_phiidCE_DC_CD = get_phiidCE_DC_CD(micro_variables, ...
		measure, method, time_lag, red_func, varargin);
	
	% Function description: get_phiidCE_DC_C() calculates PhiID-CE, 
	% PhiID-DC, and PhiID-CD for any micro variable specified 
	% in [micro_variables], and the parameters given; it will return 
	% the PhiID-measures as specified in [measure].
	
	% Inputs:	
	%
	% Required:	micro_variables			1x1 struct with micro 
	%							variable names in fields, 
	%							and micro variable 
	%							time-series as values
	%		measure				cell array
	%		method				character array
	%		time_lag				double
	%
	% Optional: kraskov_param			double
	%
	% Outputs:	all_phiidCE_DC_CD			1 x n struct with with micro
	%							variables in fields, and 
	%							PhiID-measures as values 
	%							(n refers to number of PhiID-
	%							measures specified in [measure])
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'micro_variables', @isstruct);
	addRequired(p,'measure', @iscell);
	addRequired(p,'method', @ischar);
	addRequired(p,'time_lag', @isdouble);
	addRequired(p,'red_func', @ischar);
	
	% optional name-value pair variables: 
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, @isdouble);
	
	default_bin_number = 3;
	addParameter(p,'bin_number', default_bin_number, @isdouble);
	
	parse(p, micro_variables, measure, method, time_lag, ...
		red_func, varargin{:});
	
	measure				= p.Results.measure;
	method				= p.Results.method;
	time_lag				= p.Results.time_lag;
	red_func				= p.Results.red_func;
	micro_variables			= p.Results.micro_variables;
	kraskov_param			= p.Results.kraskov_param;
	bin_number				= p.Results.bin_number;

	% get number of micro variables
	n_micro_variables = length(fieldnames(micro_variables)); 
	
	% get names of micro variables
	fieldnames_micro = fieldnames(micro_variables);
	
	% loop over all micro and macro variables and calculate practical CE
	for i = 1:n_micro_variables;
		
		micro = micro_variables.(fieldnames_micro{i});
			
		phiid_atoms = get_phiid_atoms(micro, method, time_lag, red_func);
		
		phiidCE = phiid_atoms.str + phiid_atoms.stx + ...
			phiid_atoms.sty + phiid_atoms.sts;
		
		phiidDC = phiid_atoms.str + phiid_atoms.stx + phiid_atoms.sty;
		
		phiidCD = phiidCE - phiidDC;
		
		% calculate PhiID-CE, PhiID-DC, and PhiID-CE for a given micro-macro
		% combination, and store results in struct
		
		all_phiidCE.(fieldnames_micro{i}) = phiidCE;
		all_phiidDC.(fieldnames_micro{i}) = phiidDC;
		all_phiidCD.(fieldnames_micro{i}) = phiidCD;
		
	end 
	
	% put all structs of different PhiID-measures into one array, in the order as 
	% specified in [measure]
	
	all_phiidCE_DC_CD = [];
	for p = 1:length(measure)
		
		if any(strcmp(measure, 'phiidCE'))
			all_phiidCE_DC_CD.all_phiidCE = all_phiidCE;
		end
		if any(strcmp(measure, 'phiidDC'))
			all_phiidCE_DC_CD.all_phiidDC = all_phiidDC;
		end
		if any(strcmp(measure, 'phiidCD'))
			all_phiidCE_DC_CD.all_phiidCD = all_phiidCD;
		end 
	end 
end 
