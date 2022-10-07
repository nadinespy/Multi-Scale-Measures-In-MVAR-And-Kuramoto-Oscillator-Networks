function phiid_atoms = get_phiid_atoms(micro, method, time_lag, red_func, varargin)
	
	% Function description: get_phiid_atoms() calculates PhiID, 
	% PhiID-DC, and PhiID-CD for any micro variable specified 
	% in [micro_variables], and the parameters given; it will return 
	% the PhiID-measures as specified in [measure].
	
	% Inputs:	
	%
	% Required:	micro					double of size 
	%							[vars in micro x time_length] 
	%							time-series as values
	%		method				character array
	%		time_lag				double
	%		red_func				character array
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
	addRequired(p,'micro', @isdouble);
	addRequired(p,'method', @ischar);
	addRequired(p,'time_lag', @isdouble);
	addRequired(p,'red_func', @ischar);
	
	% optional name-value pair variables: 
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, @isdouble);
	
	parse(p, micro, method, time_lag, red_func, varargin{:});

	method				= p.Results.method;
	time_lag				= p.Results.time_lag;
	red_func				= p.Results.red_func;
	micro					= p.Results.micro;
	kraskov_param			= p.Results.kraskov_param;
	
	if strcmp(lower(method), 'gaussian');
		phiid_atoms = PhiIDFullContinuous(micro, time_lag, red_func, method);
	elseif strcmp(lower(method), 'discrete');
		phiid_atoms = PhiIDFullDiscrete(micro, time_lag, red_func);
	elseif strcmp(lower(method), 'kraskov');
		phiid_atoms = PhiIDFullContinuous(micro, time_lag, red_func, method, kraskov_param);
	end 
	
end 