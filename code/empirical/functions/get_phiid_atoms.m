function phiid_atoms = get_phiid_atoms(micro, method, time_lag, red_func, varargin)
% get_phiid_atoms() calculates PhiID for a micro variable specified in [micro].
%
% Takes as inputs a method ([method]), a redundancy function ([red_func]), and
% a time-lag ([time_lag]), and will return all 16 PhiID-atoms in a structure.
%	
% Example: phiid_atoms = get_phiid_atoms(micro, method, time_lag, red_func, ...
%	     'kraskov_param', kraskov_param)
%
% INPUTS - required:	
%    micro_variables	    -			1x1 struct with micro 
%    micro			    -             double of size 
%							[vars in micro x time_length] 
%							time-series as values
%    method			    -			character array
%    time_lag		    -			double
%    red_func		    -			character array
%
% INPUTS - optional: 
%    kraskov_param	    -			double
%
% OUTPUT:	
%    phiid_atoms		    -			struct with 16 fields corresponding
%							to atoms 
	
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

	method		= p.Results.method;
	time_lag		= p.Results.time_lag;
	red_func		= p.Results.red_func;
	micro			= p.Results.micro;
	kraskov_param	= p.Results.kraskov_param;
	
	if strcmp(lower(method), 'gaussian');
		phiid_atoms = PhiIDFullContinuous(micro, time_lag, red_func, method);
	elseif strcmp(lower(method), 'discrete');
		phiid_atoms = PhiIDFullDiscrete(micro, time_lag, red_func);
	elseif strcmp(lower(method), 'kraskov');
		phiid_atoms = PhiIDFullContinuous(micro, time_lag, red_func, method, kraskov_param);
	end 
	
end 