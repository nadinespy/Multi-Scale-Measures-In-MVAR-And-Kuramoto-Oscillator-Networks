%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply number of variables n, residuals correlation coefficient r
% and random seed s before running this script; e.g.:
%
% >> n = 10; r = 0.8; non_monotonic_integration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('seed','var'), seed = 'shuffle'; end
rng(seed) % seed random number generator


%% calculate integration across range of NOISE CORRELATIONS and COUPLINGS, FIXED C-FACTOR

% {
Iv = zeros(length(noise_corrs),length(two_node_couplings)); % process total variance
Ig = zeros(length(noise_corrs),length(two_node_couplings)); % process generalised variance (includes covariances)
I = zeros(length(noise_corrs),length(two_node_couplings));  % integration (process multi-information)
J = zeros(length(noise_corrs),length(two_node_couplings));  % residuals multi-information
for j = 1:length(noise_corrs)
	% create a uniform correlation matrix with the same bivariate
	% correlation between every pair of variables
	if strcmp(signs_noise_corrs, 'positive')
		% positive correlation
		R = noise_corrs(j)*ones(n);
	elseif strcmp(signs_noise_corrs, 'negative')
		% negative correlation
		R = -noise_corrs(j)*ones(n);
	elseif strcmp(signs_noise_corrs, 'mixed')
		% mixed correlations
		coefficients = [-1, 1];
		r = randi([1, 2], 1);
		temp_coeff = coefficients(r);
		R = temp_coeff*noise_corrs(j)*ones(n);
	end
	R(1:n+1:n^2) = 1; % set diagonal entries to 1
	
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
