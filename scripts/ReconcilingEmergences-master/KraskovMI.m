function [ mi ] = KraskovMI(X, Y, kraskov_param)

%	Function description: KraskovMI() calculates mutual information (MI) 
%	between two continuous, (non-Gaussian) variables X and Y, or MI 
%     between one continuous (X or Y) and one discrete variable (X or Y). 
%
%	Inputs:	
%
%	Required:	X				1D double array
%			Y				1D double array
%			kraskov_param		number of k-nearest neighbours
%							(int)
%
%	Outputs:	mi				mutual information between X & Y
%
% Nadine Spychala, Sep 2022 

	% use inputParser to declare variables
	p = inputParser;
	
	% required arguments:
	addRequired(p,'X', @isdouble);
	addRequired(p,'Y', @isdouble);
	addRequired(p, 'kraskov_param', @isdouble);
	
	parse(p, X, Y, kraskov_param);
	
	X					= p.Results.X;
	Y					= p.Results.Y;
	kraskov_param			= p.Results.kraskov_param;

	% parameter checks
	if ~(isvector(X) && isvector(Y) && length(X) == length(Y))
		error("X and Y must be vectors of the same length.");
	end

	%compute MI
	if isfloat(X) && isinteger(Y)
		mi = mi_discrete_cont(X, Y, kraskov_param);
	elseif isfloat(Y) && isinteger(X)
		mi = mi_discrete_cont(Y, X, kraskov_param);
	else
		mi = mi_cont_cont(X, Y, kraskov_param);
	end 

end 
