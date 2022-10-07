function [ delta, v_mi, x_mi ] = EmergenceDelta(X, V, tau, varargin)
%% EMERGENCEDELTA Compute downward causation criterion from data
%
%     DELTA = EMERGENCEDELTA(X, V) computes the downward causation criterion
%     delta for the system with micro time series X and macro time series Y.
%     Micro data must be of size TxD and macro data of size Tx1. Note that for
%     the theory to hold V has to be a (possibly stochastic) function of X.
%
%     DELTA = EMERGENCEDELTA(X, V, TAU) uses a time delay of TAU samples to
%     compute time-delayed mutual information. (default: 1)
%
%     DELTA = EMERGENCEDELTA(X, V, TAU, METHOD) estimates mutual info using a
%     particular METHOD. Can be 'discrete' or 'gaussian' or 'kraskov'. 
%	If empty, it will try to infer the most suitable method.
%
%     DELTA = EMERGENCEDELTA(X, V, TAU, METHOD, KRASKOV_PARAM) estimates 
%     mutual info using KRASKOV_PARAM. 
%
%     [DELTA, V_MI, X_MI] = EMERGENCEDELTA(X, V, ...) also returns the mutual
%     info in macro and micro variables, such that DELTA = max(V_MI - X_MI).
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
	
	% parameter checks and initialisation
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
	
	
	% compute mutual infos and delta
	
	if strcmp(lower(method), 'kraskov')
		
		v_mi = arrayfun(@(j) MI_fun(V(1:end-tau), X(1+tau:end,j), ...
			kraskov_param), 1:size(X,2));
		x_mi = arrayfun(@(j) sum(arrayfun(@(i) MI_fun(X(1:end-tau,i), ...
			X(1+tau:end,j), kraskov_param), ...
			1:size(X,2))), 1:size(X,2));
		
	else
		
		v_mi = arrayfun(@(j) MI_fun(V(1:end-tau), X(1+tau:end,j)), ...
			1:size(X,2));
		x_mi = arrayfun(@(j) sum(arrayfun(@(i) MI_fun(X(1:end-tau,i), ...
			X(1+tau:end,j)), 1:size(X,2))), 1:size(X,2));
		
	end
	
	delta = max(v_mi - x_mi);
	
	
	if nargout < 2
		clearvars x_mi, v_mi;
	end

