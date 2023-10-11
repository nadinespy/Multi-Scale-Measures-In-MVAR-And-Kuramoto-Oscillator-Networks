function [ redred, localred ] = DoubleRedundancyMMIDiscrete(varargin)
%%DOUBLEREDUNDANCYMMI Compute the PhiID double-redundancy of discrete input 
% data, using the MMI PID. It uses a plug-in MI estimator to find
% the MIB.
%
% NOTE: assumes JIDT has been already added to the javaclasspath.
%
%   R = DOUBLEREDUNDANCYMMI(X, tau), where X is a D-by-T data matrix of D
%   dimensions for T timesteps, and TAU is an integer integration timescale,
%   computes the double-redundancy across the minimum information bipartition
%   (MIB) of X. If TAU is not provided, it is set to 1.
%
%   R = DOUBLEREDUNDANCYMMI(X, tau, estimator), where X is a D-by-T 
%   data matrix of D dimensions for T timesteps, and TAU is an integer integration 
%   timescale, ESTIMATOR is a character array denoting the estimator used to 
%   estimate mutual information. Implemented estimators are 'PlugIn' and 'NSB' 
%   (default: 'NSB'), computes the double-redundancy across the minimum 
%   information bipartition (MIB) of X. If TAU is not provided, it is set to 1.
%
%   R = DOUBLEREDUNDANCYMMI('X1', X1, 'X2', X2, 'Y1', Y1, 'Y2', Y2), where 
%   all inputs are matrices with the same number of columns (i.e. same number 
%   of samples), computes the double-redundancy of the mutual info between 
%   them, I(X1, X2; Y1, Y2).
%
%   [R, L] = DOUBLEREDUNDANCYMMI(...) returns the local double-redundancy
%   values for each sample in the input. (NOTE: not available for the NSB 
%   estimator.
%
% If input data is discrete-compatible (as per ISDISCRETE), it is passed
% directly to the underlying information-theoretic calculators. If it isn't
% (e.g. if it is real-valued data), it is mean-binarised first.
%
% Reference:
%   Mediano*, Rosas*, Carhart-Harris, Seth and Barrett (2019). Beyond
%   Integrated Information: A Taxonomy of Information Dynamics Phenomena.
%
% Pedro Mediano, Jan 2021
% adapted by Nadine Spychala, May 2023

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
	red_func    = p.Results.red_func;
	
	if nargin <= 4 && (isempty(X) == false)
		atoms = private_TDMMI(X, 'tau', tau, 'red_func', red_func);
		
	elseif nargin >= 8 && (isempty(X) == true)
		atoms = private_FourVectorMMI(X1, X2, Y1, Y2, 'red_func', red_func);
	else
		error('Wrong number of arguments. See `help DoubleRedundancyMMI` for help.');
	end
	
	redred = mean(atoms(isfinite(atoms)));
	
	if nargout > 1
		localred = atoms;
	end

end


%*********************************************************
%*********************************************************
function [ redred ] = private_TDMMI(X, varargin)

	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required arguments:
	addRequired(p,'X', @isdouble);

	% optional name-value variables:
	default_tau = 1;
	addParameter(p,'tau', @isdouble);
	
	default_red_func = 'mmi';
	addParameter(p,'red_func', default_red_func, @ischar);
	
	parse(p, varargin{:});
	
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

	integer_tau = ~isinf(tau) & floor(tau) == tau;
	if ~integer_tau || tau < 1
		error('Timescale tau needs to be a positive integer.');
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
	
	redred = private_FourVectorMMI(X(p1, 1:end-tau), X(p2, 1:end-tau), ...
		X(p1, 1+tau:end), X(p2, 1+tau:end), 'red_func', red_func);

end


%*********************************************************
%*********************************************************
function [ redred ] = private_FourVectorMMI(X1, X2, Y1, Y2, varargin)
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required arguments
	addRequired(p,'X1', @isdouble);
	addRequired(p,'X2', @isdouble);
	addRequired(p,'Y1', @isdouble);
	addRequired(p,'Y2', @isdouble);
	
	% optional name value pairs:
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
	binarify = @(v) isdiscrete(v)*v + (~isdiscrete(v))*(v > mean(v, 2));
	src = {binarify(X1), binarify(X2)};
	tgt = {binarify(Y1), binarify(Y2)};
	
	% take double-redundancy as the minimum MI between either src or tgt
	redred = inf([T, 1]);
	for i=1:length(src)
		for j=1:length(tgt)
			
			% Bayes estimation doesn't work
			if strcmp(red_func, 'bmmi')
				mi = QuasiBayesMI(src{i}, tgt{j});
			else
				x = ensure_combined(src{i});
				y = ensure_combined(tgt{j});
				
				miCalc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', max(x)+1, max(y)+1, 0);
				miCalc.initialise();
				miCalc.addObservations(x', y');
				mi = miCalc.computeLocalFromPreviousObservations(x', y');
			end
			
			if mean(mi) < mean(redred)
				redred = mi;
			end
		end
	end

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


