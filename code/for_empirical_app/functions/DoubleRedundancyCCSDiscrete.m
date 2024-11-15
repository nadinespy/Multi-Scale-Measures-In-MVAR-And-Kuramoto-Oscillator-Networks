function [ redred, localred ] = DoubleRedundancyCCSDiscrete(varargin)
%%DOUBLEREDUNDANCYCCS Compute the PhiID double-redundancy of discrete 
% input data using the CCS PID. It uses a plug-in MI estimator.
%
% NOTE: assumes JIDT has been already added to the javaclasspath.
%
%   R = DOUBLEREDUNDANCYCCS(X, tau), where X is a D-by-T data matrix of D
%   dimensions for T timesteps, and TAU is an integer integration timescale,
%   computes the double-redundancy across the minimum information bipartition
%   (MIB) of X. If TAU is not provided, it is set to 1.
%
%   R = DOUBLEREDUNDANCYCCS(X, tau, method, kraskov_param), where X is a D-by-T 
%   data matrix of D dimensions for T timesteps, and TAU is an integer integration 
%   timescale, METHOD is a character array denoting the method to use (Kraskov or 
%   Gaussian), KRASKOVPARAM is an integer indicating k-nearest neighbours,
%   computes the double-redundancy across the minimum information bipartition
%   (MIB) of X. If TAU is not provided, it is set to 1.
%
%   R = DOUBLEREDUNDANCYCCS('X1', X1, 'X2', X2, 'Y1', Y1, 'Y2', Y2), where 
%   all inputs are matrices with the same number of columns (i.e. same number 
%   of samples), computes the double-redundancy of the mutual info between 
%   them, I(X1, X2; Y1, Y2).
%
%   [R, L] = DOUBLEREDUNDANCYCCS(...) returns the local double-redundancy
%   values for each sample in the input.
% 
% If input data is discrete-compatible (as per ISDISCRETE), it is passed
% directly to the underlying information-theoretic calculators. If it isn't
% (e.g. if it is real-valued data), it is mean-binarised first.
%
% Reference:
%   Mediano*, Rosas*, Carhart-Harris, Seth and Barrett (2019). Beyond
%   Integrated Information: A Taxonomy of Information Dynamics Phenomena.
%
%   Ince (2016). Measuring multivariate redundant information with pointwise
%   common change in surprisal
%
% Pedro Mediano, Jan 2021
% Adapted by Nadine Spychala, May 2023

	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% optional name-value pair variables:
	default_X = [];
	addParameter(p,'X', default_X, @isdouble);
	
	default_tau = 1;
	addParameter(p,'tau', default_tau, @isdouble);

	default_red_func = 'mmi';
	addParameter(p,'red_func', default_red_func, @ischar);
	
	default_X1 = [];
	addParameter(p,'X1', default_X1, @isdouble);
	default_X2 = [];
	addParameter(p,'X2', default_X2, @isdouble);
	default_Y1 = [];
	addParameter(p,'Y1', default_Y1, @isdouble);
	default_Y2 = [];
	addParameter(p,'Y2', default_Y2, @isdouble);
	
	parse(p, varargin{:});
	
	X		= p.Results.X;
	X1		= p.Results.X1;
	X2		= p.Results.X2;
	Y1		= p.Results.Y1;
	Y2		= p.Results.Y2;
	tau		= p.Results.tau;
	red_func     = p.Results.red_func;
		
	if nargin <= 4 && (isempty(X) == false)
		R = private_TDCCS(X, 'tau', tau, 'red_func', red_func);
		
	elseif nargin >= 8 && (isempty(X) == true)
		R = private_FourVectorCCS(X1, X2, Y1, Y2, 'red_func', red_func);
	else
		error('Wrong number of arguments. See `help DoubleRedundancyMMI` for help.');
	end
	
	redred = mean(R(isfinite(R)));
	
	if nargout > 1
		localred = R;
	end

end


%*********************************************************
%*********************************************************
function [ redred ] = private_TDCCS(X, varargin)
		  
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required arguments:
	addRequired(p,'X', @isdouble);
	
	% optional name-value variables:
	default_tau = 1;
	addParameter(p, 'tau', default_tau, @isdouble);
	
	default_red_func = 'mmi';
	addParameter(p,'red_func', default_red_func, @ischar);

	parse(p, X, varargin{:});
	
	X		= p.Results.X;
	tau		= p.Results.tau;
	red_func    = p.Results.red_func;

	% argument checks and parameter initialisation
	if isempty(X) || ~ismatrix(X)
	  error('Input must be a 2D data matrix');
	end
	
	[D, T] = size(X);
	if T <= D
	  error(sprintf(['Your matrix has %i dimensions and %i timesteps. ', ...
		  'If this is true, you cant compute a reasonable covariance matrix. ', ...
		  'If it is not true, you may have forgotten to transpose the matrix'], D, T));
	end
	
	% binarise the data, if not already discrete
	if ~isdiscrete(X)
		X = 1*(X > mean(X, 2));
	end
	
	% use JIDT to compute Phi and MIB
	phiCalc = javaObject('infodynamics.measures.discrete.IntegratedInformationCalculatorDiscrete', 2, size(X, 1));
	if tau > 1
		phiCalc.setProperty(phiCalc.PROP_TAU, num2str(tau));
	end
	phiCalc.setObservations(octaveToJavaIntMatrix(X'));
	phi = phiCalc.computeAverageLocalOfObservations();
	mib = phiCalc.getMinimumInformationPartition();
	
	% extract MIB partition indices
	p1 = str2num(mib.get(0).toString()) + 1;
	p2 = str2num(mib.get(1).toString()) + 1;
	
	% here, we define X1, X2, Y1, and Y2, so the partitions at time t (X1, X2),
	% and the partitions at time t+1 (Y1, Y2)
	redred = private_FourVectorCCS(X(p1,1:end-tau), X(p2,1:end-tau), ...
			X(p1,1+tau:end), X(p2,1+tau:end), 'red_func', red_func);

end


%*********************************************************
%*********************************************************
function [ redred ] = private_FourVectorCCS(X1, X2, Y1, Y2, varargin)

	% use inputParser to declare variables
	p = inputParser;
	
	% required arguments:
	addRequired(p,'X1', @isdouble);
	addRequired(p,'X2', @isdouble);
	addRequired(p,'Y1', @isdouble);
	addRequired(p,'Y2', @isdouble);

	% optional name-value variables
	default_red_func = 'mmi';
	addParameter(p,'red_func', default_red_func, @ischar);

	parse(p, X1, X2, Y1, Y2, varargin{:});
	
	X1		= p.Results.X1;
	X2		= p.Results.X2;
	Y1		= p.Results.Y1;
	Y2		= p.Results.Y2;
	red_func    = p.Results.red_func;
	
	% argument checks and parameter initialisation
	T = size(X1, 2);
	if size(X2, 2) ~= T || size(Y1, 2) ~= T || size(Y2, 2) ~= T
		error('All input vectors must have the same length');
	end
	
	% binarise data (if not already discrete) and stack for easier handling
	binarify = @(v) ensure_combined(isdiscrete(v)*v + (~isdiscrete(v))*(v > mean(v, 2)));
	X = [binarify(X1); binarify(X2); binarify(Y1); binarify(Y2)];

	Ixytab = localmi(X, [1, 2], [3, 4]);
	
	Ixta = localmi(X, 1, 3);
	Ixtb = localmi(X, 1, 4);
	Iyta = localmi(X, 2, 3);
	Iytb = localmi(X, 2, 4);
	
	Ixyta = localmi(X, [1, 2], 3);
	Ixytb = localmi(X, [1, 2], 4);
	Ixtab = localmi(X, 1, [3, 4]);
	Iytab = localmi(X, 2, [3, 4]);
	
	Rxytab = localred(Ixtab, Iytab, Ixytab);
	Rabtxy = localred(Ixyta, Ixytb, Ixytab);
	Rxyta  = localred(Ixta, Iyta, Ixyta);
	Rxytb  = localred(Ixtb, Iytb, Ixytb);
	Rabtx  = localred(Ixta, Ixtb, Ixtab);
	Rabty  = localred(Iyta, Iytb, Iytab);
	
	% this quantity equals redred - synsyn
	double_coinfo = - Ixta - Ixtb - Iyta - Iytb + ...
		+ Ixtab + Iytab + Ixyta + Ixytb - Ixytab + ...
		+ Rxyta + Rxytb - Rxytab + ...
		+ Rabtx + Rabty - Rabtxy;
	
	% double_coinfo above is equal to:
	%	   - I(1(t);1(t+1)) - I(1(t);2(t+1)) - I(2(t);1(t+1)) - I(2(t);2(t+1))
	%	   + I(1(t);1(t+1),2(t+1)) + I(2(t);1(t+1),2(t+1)) + I(1(t),2(t);1(t+1)) 
	%	   + I(1(t),2(t);2(t+1)) - I(1(t),1(t+1);2(t),2(t+1))
	%	   + a bunch of redundancy terms derived for a univariate target
	
	redred = all(~diff([sign(Ixta), sign(Ixtb), sign(Iyta), sign(Iytb), sign(double_coinfo)],1,2),2).*double_coinfo;

end


%*********************************************************
%*********************************************************
function [ V ] = ensure_combined(U)
	if size(U, 1) == 1
		V = U;
	else
		[~,~,V] = unique(U', 'rows');
		V = V(:)' - 1;
	end
end

function [ l ] = localmi(X, src, tgt)
	
	x = ensure_combined(X(src,:));
	y = ensure_combined(X(tgt,:));
	
	miCalc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', max(x)+1, max(y)+1, 0);
	miCalc.initialise();
	miCalc.addObservations(x', y');
	l = miCalc.computeLocalFromPreviousObservations(x', y');

end

function [ l ] = localred(mi1, mi2, mi12)
	
	% Co-information (equivalent to interaction information): defined as
	% I(S1;T) + I(S2;T) - I(S1,S2;T) for sources S1, S2, and target T.
	% Negative values denote synergy, positive ones redundancy. We flip
	% signs for positive values to denote synergy such that
	% co-information = I(S1,S2;T) - I(S1;T) - I(S2;T)
	
	c = mi12 - mi1 - mi2;
  
	% store signs of all terms (sign of c is now flipped for positive 
	% values to indicate redundancy
	
	signs = [sign(mi1), sign(mi2), sign(mi12), sign(-c)];
  
	% all() determines if the elements are all nonzero or logical 1 
	% (true); 
	% signs == signs(:,1) --> creates new matrix indicating whether 
	% columns in signs are equal to the first column of signs (all 
	% (elements of first column in the new matrix would trivially be 1);
	% all(signs == signs(:,1), 2) --> determines whether all elements of
	% a given row in the new matrix are nonzero or 1 (yields column 
	% vector with ones and zeros);
	% multiply column vector with (-c) (yields zero for cases where signs 
	% had not been equal)  
	
	l = all(signs == signs(:,1), 2).*(-c);
end

