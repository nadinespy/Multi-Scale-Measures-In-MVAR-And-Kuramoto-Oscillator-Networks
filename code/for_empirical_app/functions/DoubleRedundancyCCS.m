function [ redred, localred ] = DoubleRedundancyCCS(varargin)
%%DOUBLEREDUNDANCYCCS Compute the PhiID double-redundancy of input data,
% assuming it follows a Gaussian distribution and using the CCS PID.
%
% NOTE: assumes JIDT has been already added to the javaclasspath.
%
%   R = DOUBLEREDUNDANCYCCS('X', X, 'tau', tau), where X is a D-by-T data matrix of D
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

	default_method = 'Gaussian';
	addParameter(p,'method', default_method, @ischar);
	
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, @isdouble);
	
	default_X1 = [];
	addParameter(p,'X1', default_X1, @isdouble);
	default_X2 = [];
	addParameter(p,'X2', default_X2, @isdouble);
	default_Y1 = [];
	addParameter(p,'Y1', default_Y1, @isdouble);
	default_Y2 = [];
	addParameter(p,'Y2', default_Y2, @isdouble);
	
	parse(p, varargin{:});
	
	X			= p.Results.X;
	X1			= p.Results.X1;
	X2			= p.Results.X2;
	Y1			= p.Results.Y1;
	Y2			= p.Results.Y2;
	tau			= p.Results.tau;
	method		= p.Results.method;
	kraskov_param	= p.Results.kraskov_param;
	
	if nargin == 2
		R = private_TDCCS(X);
	elseif nargin == 4
		R = private_TDCCS(X, 'tau', tau);
	elseif nargin == 6
		R = private_TDCCS(X, 'tau', tau, 'method', method);
	elseif nargin == 6 && (isempty(X) == false)
		R = private_TDCCS(X, 'tau', tau, 'method', method, 'kraskov_param', kraskov_param);
	elseif nargin == 8 && (isempty(X) == true)
		R = private_FourVectorCCS(X1, X2, Y1, Y2);
	elseif nargin == 10 
		R = private_FourVectorCCS(X1, X2, Y1, Y2, 'method', method);
	elseif nargin == 12
		R = private_FourVectorCCS(X1, X2, Y1, Y2, 'method', method, 'kraskov_param', kraskov_param);
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

	% optional name-value pair variables:
	default_tau = 1;
	addParameter(p,'tau', @isdouble);
	
	default_method = 'Gaussian';
	addParameter(p,'method', default_method, ...
		@ischar);
	
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, ...
		@isdouble);
	
	parse(p, varargin{:});
	
	X			= p.Results.X;
	tau			= p.Results.tau;
	method		= p.Results.method;
	kraskov_param	= p.Results.kraskov_param;
	
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
	
	% Create copy of the data scaled to unit variance (for numerical stability)
	sX = X./repmat(std(X')', [1, T]);
	
	
	% use JIDT to compute Phi and MIB
	phiCalc = javaObject('infodynamics.measures.continuous.gaussian.IntegratedInformationCalculatorGaussian');
	
	if tau > 1
		phiCalc.setProperty(phiCalc.PROP_TAU, num2str(tau));
	end
	
	phi = phiCalc.compute(octaveToJavaDoubleMatrix(sX'));
	mib = phiCalc.getMinimumInformationPartition();
	
	% extract MIB partition indices
	p1 = str2num(mib.get(0).toString()) + 1;
	p2 = str2num(mib.get(1).toString()) + 1;
	
	% here, we define X1, X2, Y1, and Y2, so the partitions at time t (X1, X2),
	% and the partitions at time t+1 (Y1, Y2)
	if strcmp(lower(method), 'kraskov');
		
		redred = private_FourVectorCCS(X(p1,1:end-tau), X(p2,1:end-tau), ...
			X(p1,1+tau:end), X(p2,1+tau:end), 'method', method, 'kraskov_param', kraskov_param);
		
	else redred = private_FourVectorCCS(X(p1,1:end-tau), X(p2,1:end-tau), ...
			X(p1,1+tau:end), X(p2,1+tau:end), 'method', method);
		
	end
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
	
	% optional name-value variables:
	default_method = 'Gaussian';
	addParameter(p, 'method', default_method, @ischar);
	default_kraskov_param = 3;
	addParameter(p, 'kraskov_param', default_kraskov_param, @isdouble);

	parse(p, X1, X2, Y1, Y2, varargin{:});
	
	X1			= p.Results.X1;
	X2			= p.Results.X2;
	Y1			= p.Results.Y1;
	Y2			= p.Results.Y2;
	method		= p.Results.method;
	kraskov_param	= p.Results.kraskov_param;
	
	if strcmp(lower(method), 'kraskov');
		error('Kraskov estimation not implemented yet.');
	end
	
	% argument checks and parameter initialisation
	T = size(X1, 2);
	if size(X2, 2) ~= T || size(Y1, 2) ~= T || size(Y2, 2) ~= T
		error('All input vectors must have the same length');
	end
	
	% indices of source and target variables in the joint covariance matrix
	p1 = 1:size(X1, 1);                     % indices of first source partition (e.g., in an 8-element system [1,2,3,4])
	p2 = p1(end)+1:p1(end)+size(X2, 1);     % indices of second source partition (e.g., in an 8-element system [5,6,7,8])
	t1 = p2(end)+1:p2(end)+size(Y1, 1);     % indices of first target partition (e.g., " [9,10,11,12])
	t2 = t1(end)+1:t1(end)+size(Y2, 1);     % indices of second target partition (e.g., " [13,14,15,16])
	D = t2(end);                            % total number of indices
	
	% stack data for easier handling (also scale to unit variance for numerical stability)
	X = [X1; X2; Y1; Y2];                   % stack all variables from partitions row-wise
	sX = X./repmat(std(X')', [1, T]);       % we again scale to unit variance, because the matrix has changed
	
	% compute mean and covariance for all the data (to be used by the local IT functions below)
	% S: e.g., in an 8-element system, the size will be 16x16, it's the time-lagged covariance matrix
	S = cov(sX');                                          
	mu = mean(sX');
	assert(all(size(mu) == [1, D]) && all(size(S) == [D, D]));
	
	% define local information-theoretic functions (lacking: Kraskov version of estimating non-Gaussian multivariate entropy)
	% local multivariate entropy, h() takes as an input the indices of the variables to consider in sX
	h = @(idx) -log(mvnpdf(sX(idx,:)', mu(idx), S(idx,idx)));
	
	% sample from estimated Gaussian and do Monte Carlo averaging - don't know anymore what this was good for
	% nb_samples = 5000;
	% % mvnrnd() returns an m-by-d matrix of random vectors sampled from m separate d-dimensional normal distributions
	% x = mvnrnd(mu, S, nb_samples); 
	
	% pre-compute entropies necessary for all IT quantities
	h_p1 = h(p1);                      % entropy of partition 1 at t           H(1(t))
	h_p2 = h(p2);                      % entropy pf partiiton 2 at t           H(2(t))
	h_t1 = h(t1);                      % entropy of partition 1 at t+1         H(1(t+1))
	h_t2 = h(t2);                      % entropy of partition 2 at t+1         H(2(t+1))
	
	h_p1p2 = h([p1 p2]);               % multivariate entropy (ME) of partition 1 & 2 at t         H(1(t),   2(t))
	h_t1t2 = h([t1 t2]);               % ME of partition 1 & 2 at t+1                              H(1(t+1), 2(t+1))
	h_p1t1 = h([p1 t1]);               % ME of partition 1 at t & t+1                              H(1(t),   1(t+1))
	h_p1t2 = h([p1 t2]);               % ME of partition 1 at t & partition 2 at t+1               H(1(t),   2(t+1))
	h_p2t1 = h([p2 t1]);               % ME of partition 2 at t & partition 1 at t+1               H(2(t),   1(t+1))
	h_p2t2 = h([p2 t2]);               % ME of partition 2 at t & t+1                              H(2(t),   2(t+1))
	
	h_p1p2t1 = h([p1 p2 t1]);          % ME of partition 1 & 2 at t & partition 1 at t+1           H(1(t),   2(t),     1(t+1))
	h_p1p2t2 = h([p1 p2 t2]);          % ME of partition 1 & 2 at t & partition 2 at t+1           H(1(t),   2(t),     2(t+1))
	h_p1t1t2 = h([p1 t1 t2]);          % ME of partition 1 at t & t+1 & partition 2 at t+1         H(1(t),   1(t+1),   2(t+1))
	h_p2t1t2 = h([p2 t1 t2]);          % ME of partition 2 at t & t+1 & partition 1 at t           H(2(t),   2(t+1),   1(t+1))
	
	h_p1p2t1t2 = h([p1 p2 t1 t2]);     % ME of partition 2 at t & t+1 & partition 1 at t & t+1     H(2(t),   2(t+1),   1(t),     2(t+1))
	
	% compute PhiID quantities as entropy combinations
	
	% Ixytab: all 16 atoms (this is what we're trying to decompose)
	Ixytab = h_p1p2 + h_t1t2 - h_p1p2t1t2;          % I(1(t),1(t+1);2(t),2(t+1))    H(1(t),2(t)) + H(1(t+1),2(t+1)) - H(1(t),1(t+1),2(t),2(t+1))          
	
	Ixta = h_p1 + h_t1 - h_p1t1;                    % I(1(t);1(t+1))                H(1(t))      + H(1(t+1))        - H(1(t),1(t+1))
	Ixtb = h_p1 + h_t2 - h_p1t2;                    % I(1(t);2(t+1))                H(1(t))      + H(2(t+1))        - H(1(t),2(t+1))
	Iyta = h_p2 + h_t1 - h_p2t1;                    % I(2(t);1(t+1))                H(2(t))      + H(1(t+1))        - H(2(t),1(t+1))
	Iytb = h_p2 + h_t2 - h_p2t2;                    % I(2(t);2(t+1))                H(2(t))      + H(2(t+1))        - H(2(t),2(t+1))
	
	Ixyta = h_p1p2 + h_t1 - h_p1p2t1;               % I(1(t),2(t);1(t+1))           H(1(t),2(t)) + H(1(t+1))        - H(1(t),2(t),1(t+1))
	Ixytb = h_p1p2 + h_t2 - h_p1p2t2;               % I(1(t),2(t);2(t+1))           H(1(t),2(t)) + H(2(t+1))        - H(1(t),2(t),2(t+1))
	Ixtab = h_p1 + h_t1t2 - h_p1t1t2;               % I(1(t);1(t+1),2(t+1))         H(1(t))      + H(1(t+1),2(t+1)) - H(1(t),1(t+1),2(t+1))
	Iytab = h_p2 + h_t1t2 - h_p2t1t2;               % I(2(t);1(t+1),2(t+1))         H(2(t))      + H(1(t+1),2(t+1)) - H(2(t),1(t+1),2(t+1))
	
	Rxytab = localred(Ixtab, Iytab, Ixytab);        % I(1(t),2(t);1(t+1))          - I(2(t);1(t+1))        - I(1(t);1(t+1))
	Rabtxy = localred(Ixyta, Ixytb, Ixytab);        % I(1(t),2(t);2(t+1))          - I(1(t);2(t+1))	       - I(2(t);2(t+1))
	Rxyta  = localred(Ixta, Iyta, Ixyta);           % I(1(t),1(t+1);2(t),2(t+1))   - I(2(t);1(t+1),2(t+1)) - I(1(t);1(t+1),2(t+1))
	Rxytb  = localred(Ixtb, Iytb, Ixytb);           % I(1(t);1(t+1),2(t+1))	       - I(1(t);1(t+1))	       - I(1(t);2(t+1))
	Rabtx  = localred(Ixta, Ixtb, Ixtab);           % I(2(t);1(t+1),2(t+1))	       - I(2(t);1(t+1))	       - I(2(t);2(t+1))
	Rabty  = localred(Iyta, Iytb, Iytab);           % I(1(t),1(t+1);2(t),2(t+1))   - I(1(t),2(t);1(t+1))   - I(1(t),2(t);2(t+1))
	
	% difference to DoubleRedundancyMMI():
	% Rxyta  = RedFun(sX, p1, p2, t1, Ixta, Iyta, Ixyta);             % take minimum of I(1(t);1(t+1))           & I(2(t);1(t+1))
	% Rxytb  = RedFun(sX, p1, p2, t2, Ixtb, Iytb, Ixytb);             % take minimum of I(2(t);2(t+1))           & I(1(t);2(t+1))
	% Rxytab = RedFun(sX, p1, p2, [t1 t2], Ixtab, Iytab, Ixytab);     % take minimum of I(1(t);1(t+1),2(t+1))    & I(2(t);1(t+1),2(t+1))
	% Rabtx  = RedFun(sX, t1, t2, p1, Ixta, Ixtb, Ixtab);             % take minimum of I(1(t);1(t+1))           & I(1(t);2(t+1))
	% Rabty  = RedFun(sX, t1, t2, p2, Iyta, Iytb, Iytab);             % take minimum of I(2(t);1(t+1))           & I(2(t);2(t+1))
	% Rabtxy = RedFun(sX, t1, t2, [p1 p2], Ixyta, Ixytb, Ixytab);     % take minimum of I(1(t),2(t);1(t+1))      & I(1(t),2(t);2(t+1))
	
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
	
	signs = [sign(Ixta), sign(Ixtb), sign(Iyta), sign(Iytb), sign(double_coinfo)];
	% I(1(t);1(t+1)), I(1(t);2(t+1)), I(2(t);1(t+1)), I(2(t);2(t+1)), redred - synsyn
	redred = all(signs == signs(:,1), 2).*double_coinfo;

end

function [ l ] = localred(mi1, mi2, mi12)
  c = mi12 - mi1 - mi2;
  signs = [sign(mi1), sign(mi2), sign(mi12), sign(-c)];
  l = all(signs == signs(:,1), 2).*(-c);
end

