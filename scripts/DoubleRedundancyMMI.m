function [ redred, localred ] = DoubleRedundancyMMI(varargin)
%%DOUBLEREDUNDANCYMMI Compute the PhiID double-redundancy of input data,
% assuming it follows a Gaussian distribution and using the MMI PID.
%
% NOTE: assumes JIDT has been already added to the javaclasspath.
%
%   R = DOUBLEREDUNDANCYMMI(X, tau), where X is a D-by-T data matrix of D
%   dimensions for T timesteps, and TAU is an integer integration timescale,
%   computes the double-redundancy across the minimum information bipartition
%   (MIB) of X. If TAU is not provided, it is set to 1.
%
%   R = DOUBLEREDUNDANCYMMI(X, tau, method, kraskov_param), where X is a D-by-T 
%   data matrix of D dimensions for T timesteps, and TAU is an integer integration 
%   timescale, METHOD is a character array denoting the method to use (Kraskov or 
%   Gaussian), KRASKOVPARAM is an integer indicating k-nearest neighbours,
%   computes the double-redundancy across the minimum information bipartition
%   (MIB) of X. If TAU is not provided, it is set to 1.
%
%   R = DOUBLEREDUNDANCYMMI('X1', X1, 'X2', X2, 'Y1', Y1, 'Y2', Y2), where 
%   all inputs are matrices with the same number of columns (i.e. same number 
%   of samples), computes the double-redundancy of the mutual info between 
%   them, I(X1, X2; Y1, Y2).
%
%   [R, L] = DOUBLEREDUNDANCYMMI(...) returns the local double-redundancy
%   values for each sample in the input.
%
% Reference:
%   Mediano*, Rosas*, Carhart-Harris, Seth and Barrett (2019). Beyond
%   Integrated Information: A Taxonomy of Information Dynamics Phenomena.
%
% Pedro Mediano, Jan 2021
% adapted by Nadine Spychala, Oct 2022

%%
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
	
	X					= p.Results.X;
	X1					= p.Results.X1;
	X2					= p.Results.X2;
	Y1					= p.Results.Y1;
	Y2					= p.Results.Y2;
	tau					= p.Results.tau;
	method				= p.Results.method;
	kraskov_param			= p.Results.kraskov_param;

	
	if nargin == 2
		R = private_TDMMI(X);
	elseif nargin == 4
		R = private_TDMMI(X, tau);
	elseif nargin == 6
		R = private_TDMMI(X, tau, method);
	elseif nargin == 6 && (isempty(X) == false)
		R = private_TDMMI(X, tau, method, kraskov_param);
	elseif nargin == 8 && (isempty(X) == true)
		R = private_FourVectorMMI(X1, X2, Y1, Y2);
	elseif nargin == 10 
		R = private_FourVectorMMI(X1, X2, Y1, Y2, method);
	elseif nargin == 12
		R = private_FourVectorMMI(X1, X2, Y1, Y2, method, kraskov_param);
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
function [ redred ] = private_TDMMI(X, tau, varargin)

	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% optional positional arguments:
	addRequired(p,'X', @isdouble);
	addRequired(p,'tau', @isdouble);
	
	% optional name-value pair variables: 
	default_method = 'Gaussian';
	addParameter(p,'method', default_method, ...
		@ischar);
	
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, ...
		@isdouble);
	
	parse(p, varargin{:});
	
	X					= p.Results.X;
	tau					= p.Results.tau;
	method				= p.Results.method;
	kraskov_param			= p.Results.kraskov_param;
	
	
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
		
		redred = private_FourVectorMMI(X(p1,1:end-tau), X(p2,1:end-tau), ...
			X(p1,1+tau:end), X(p2,1+tau:end), method, kraskov_param);
		
	else redred = private_FourVectorMMI(X(p1,1:end-tau), X(p2,1:end-tau), ...
			X(p1,1+tau:end), X(p2,1+tau:end));
		
	end

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
	
	% optional positional arguments:
	default_method = 'Gaussian';
	addOptional(p,'method', default_method, ...
		@ischar);
	
	default_kraskov_param = 3;
	addOptional(p,'kraskov_param', default_kraskov_param, ...
		@isdouble);
	
	parse(p, X1, X2, Y1, Y2, varargin{:});
	
	X1					= p.Results.X1;
	X2					= p.Results.X2;
	Y1					= p.Results.Y1;
	Y2					= p.Results.Y2;
	method				= p.Results.method;
	kraskov_param			= p.Results.kraskov_param;
	
	% argument checks and parameter initialisation
	T = size(X1, 2);
	
	if size(X2, 2) ~= T || size(Y1, 2) ~= T || size(Y2, 2) ~= T
		error('All input vectors must have the same length');
	end
	
	
	% stack data for easier handling (also scale to unit variance for numerical stability)
	renorm = @(X) X./repmat(std(X')', [1, T]);
	src = {renorm(X1), renorm(X2)};
	tgt = {renorm(Y1), renorm(Y2)};
	
	
	% take double-redundancy as the minimum MI between either src or tgt
	redred = 1e9*ones([T, 1]);  % set to a large, but finite value
	miCalc = javaObject('infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian');
	
	for i=1:length(src)
		
		for j=1:length(tgt)
			
			s = src{i};
			t = tgt{j};	
			
			if strcmp(lower(method), 'gaussian');
				
				miCalc.initialise(size(s, 1), size(t, 1));
				miCalc.setObservations(octaveToJavaDoubleMatrix(s'), octaveToJavaDoubleMatrix(t'));
				mi = miCalc.computeLocalOfPreviousObservations();
				
				if mean(mi(isfinite(mi))) < mean(redred(isfinite(mi)))
					redred = mi;
				end
				
			elseif strcmp(lower(method), 'kraskov');
				
				mi = mi_cont_cont(s, t, kraskov_param);
				
				if mi < mean(redred(isfinite(mi)))
					redred = mi;
				end
			end
		end
	end
end


% In the loop above, we (locally) take the min between 
% I(X1(t),X2(t);X1(t+1),X2(t+1),
% I(X3(t),X4(t);X1(t+1),X2(t+1),
% I(X1(t),X2(t);X3(t+1),X4(t+1),
% I(X3(t),X4(t);X3(t+1),X4(t+1),

% if for system X, [X1, X2] give one partition, and [X3, X4] give the other.



