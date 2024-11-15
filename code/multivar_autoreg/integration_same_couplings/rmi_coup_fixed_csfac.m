%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Need MVGC2

defvar('n',				8    );	% number of variables
defvar('p',				1    );	% time-lags
defvar('w',				1    );	% decay weighting parameter (default is exponential decay of how much past terms are weighted
defvar('rmi',			1    );	% residuals' generalized correlation (multi-information)
defvar('noauto',			false);	% autocorrelation - if 'false', diagonal terms in connectivity matrix are set to zero
defvar('seed',			0    );	% random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed) % seed random number generator


%% calculate integration across range of RMI and COUPLINGS, fixed C-FACTOR

% {
Iv = zeros(length(all_rmi),length(two_node_couplings)); % process total variance
Ig = zeros(length(all_rmi),length(two_node_couplings)); % process generalised variance (includes covariances)
I = zeros(length(all_rmi),length(two_node_couplings));  % integration (process multi-information)
for j = 1:length(all_rmi);
	% create a uniform correlation matrix with the same bivariate
	% correlation between every pair of variables
	if strcmp(signs_noise_corrs, 'positive')
		% positive noise correlation
		rho = sqrt(1-exp(-all_rmi(j)));
		R = [1 rho; rho 1];
	elseif strcmp(signs_noise_corrs, 'negative')
		% negative noise correlation
		rho = sqrt(1-exp(-all_rmi(j)));
		R = [1 -rho; -rho 1];
	elseif strcmp(signs_noise_corrs, 'mixed')
		% mixed signs
		R = corr_rand(n,all_rmi(j));
	end

	% create coupling matrix by taking the ith n x n coupling matrix from coupling_matrices, 
	% and concatenate it as many times as specified by time_lags in third dimension 
	% (to account for time-lag)
	for i = 1:length(two_node_couplings)
		A = [];
		for k = 1:time_lag
			A = cat(3,A,coupling_matrices(:,:,i));
		end
		
		A = specnorm(A,1); % normalize via specnorm() - does take into account decay factor
		
		C		= var_to_autocov(csfac(csfac_index)*A,R,0);	% process covariance matrix (solve DLYAP)
		Iv(j,i)	= sum(log(diag(C)));					% process total variance
		Ig(j,i)	= log(det(C));						% process generalised variance (includes covariances)
		I(j,i)	= Iv(j,i) - Ig(j,i);					% process multi-information (nats)
	end
end
%}

