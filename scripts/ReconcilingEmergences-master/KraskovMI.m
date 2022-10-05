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
% adapted from Pedro Mediano and Fernando Rosas, Aug 2020

%%	parameter checks
	if ~(isvector(X) && isvector(Y) && length(X) == length(Y))
		error("X and Y must be vectors of the same length.");
	end

%%	compute MI

	if isfloat(X) && isdiscrete(Y)
		mi = mi_discrete_cont(X, Y, k);
	elseif isfloat(Y) && isdiscrete(X)
		mi = mi_discrete_cont(Y, X, k);
	else
		mi = mi_cont_cont(X, Y, kraskov_param);
	end 

end 
