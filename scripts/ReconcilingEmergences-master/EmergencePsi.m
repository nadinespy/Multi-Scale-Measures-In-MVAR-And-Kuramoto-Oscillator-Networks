function [ psi, v_mi, x_mi ] = EmergencePsi_2(X, V, tau, varargin)
%% EMERGENCEPSI Compute causal emergence criterion from data
%
%     PSI = EMERGENCEPSI(X, V) computes the causal emergence criterion psi for
%     the system with micro time series X and macro time series Y. Micro data
%     must be of size TxD and macro data of size Tx1. Note that for the theory
%     to hold V has to be a (possibly stochastic) function of X.
%
%     PSI = EMERGENCEPSI(X, V, TAU) uses a time delay of TAU samples to
%     compute time-delayed mutual information. (default: 1)
%
%     PSI = EMERGENCEPSI(X, V, TAU, METHOD) estimates mutual info using a
%     particular METHOD. Can be 'discrete' or 'gaussian'. If empty, it will
%     try to infer the most suitable method.
%
%     PSI = EMERGENCEPSI(X, V, TAU, METHOD, KRASKOV_PARAM) estimates mutual 
%     info using KRASKOV_PARAM.
%
%     [PSI, V_MI, X_MI] = EMERGENCEPSI(X, V, ...) also returns the mutual
%     info in macro and micro variables, such that PSI = V_MI - X_MI.
%
% Reference:
%     Rosas*, Mediano*, et al. (2020). Reconciling emergences: An
%     information-theoretic approach to identify causal emergence in
%     multivariate data. https://arxiv.org/abs/2004.08220
%
% Pedro Mediano and Fernando Rosas, Aug 2020
% Adapted by Nadine Spychala, Oct 2022

%%

	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'X', @isdouble);
	addRequired(p,'V', @isdouble);
	addRequired(p,'tau', @isdouble);
	
	% optional positional arguments:
	default_method = 'Gaussian';
	addOptional(p,'method', default_method, @ischar);
	
	% optional name-value pair variables: 
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, ...
		@isdouble);
	
	parse(p, X, V, tau, varargin{:});
	
	X					= p.Results.X;
	V					= p.Results.V;
	tau					= p.Results.tau;
	method				= p.Results.method;
	kraskov_param			= p.Results.kraskov_param;

	
	% parameter checks
	if ~isvector(V) || ~ismatrix(X)
		error("X has to be a 2D matrix and V a 1D vector.");
	end
	if length(V) ~= size(X,1)
		error("X and V must have the same length.");
	end

	if nargin < 4 || isempty(method)
		if exist('OCTAVE_VERSION', 'builtin')
			isdiscrete =  (sum(abs(X(:) - round(X(:)))) + sum(abs(V - round(V)))) < 1e-10;
		else
			isdiscrete =  iscategorical(X) || (sum(abs(X(:) - round(X(:)))) + sum(abs(V - round(V)))) < 1e-10;
		end
		if isdiscrete
			method = 'discrete';
		else
			method = 'gaussian';
		end
	end
	
	if strcmp(lower(method), 'gaussian')
		MI_fun = @GaussianMI;
	elseif strcmp(lower(method), 'discrete')
		MI_fun = @DiscreteMI;
	elseif strcmp(lower(method), 'kraskov')
		MI_fun = @KraskovMI;
	else
		error("Unknown method. Implemented options are 'gaussian' and 'discrete'.");
	end
	
	%% Compute mutual infos and psi
	
	if strcmp(lower(method), 'kraskov')
		v_mi = MI_fun(V(1:end-tau), V(1+tau:end), kraskov_param);
		x_mi = sum(arrayfun(@(j) MI_fun(X(1:end-tau,j), V(1+tau:end), kraskov_param), 1:size(X,2)));
	else
		v_mi = MI_fun(V(1:end-tau), V(1+tau:end));
		x_mi = sum(arrayfun(@(j) MI_fun(X(1:end-tau,j), V(1+tau:end)), 1:size(X,2)));
	end 

	psi = v_mi - x_mi;
			
	if nargout < 2
		clearvars x_mi v_mi;
	end
end

