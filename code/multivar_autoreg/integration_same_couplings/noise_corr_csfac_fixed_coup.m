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


%% calculate integration across range of NOISE CORRELATIONS and C-FACTORS, FIXED COUPLING

% choose the n x n coupling matrix with index coup_index from coupling_matrices, 
% and concatenate it as many times as specified by time_lags in third dimension 
% (to account for time-lag)
A = [];
for i = 1:time_lag+1
	A = cat(3,A,coupling_matrices(:,:,coup_index));
end	
A = specnorm(A,1); % specnorm() does take into account decay factor (is it correct that greater time-lag is weighed more?)

Iv = zeros(length(noise_corrs),length(csfac)); % process total variance
Ig = zeros(length(noise_corrs),length(csfac)); % process generalised variance (includes covariances)
I = zeros(length(noise_corrs),length(csfac));  % integration (process multi-information)
J = zeros(length(noise_corrs),length(csfac));  % residuals multi-information
for j = 1:length(noise_corrs);
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
	
	for i = 1:length(csfac)
		% C = dlyap(csfac(i)*A,R);                 % process covariance matrix (solve DLYAP)
		C = var_to_autocov(csfac(i)*A,R,time_lag);
		I(j,i) = sum(log2(diag(C(:,:,1)))) - log2(det(C(:,:,1))); % process multi-information (bits)
		
		% PhiID MMI & CCS
		full_time_lagged_cov = [C(:,:,1), C(:,:,time_lag+1); C(:,:,time_lag+1), C(:,:,1)];
		phiid_atoms_MMI = PhiIDFull_Analytical(full_time_lagged_cov, 'MMI');
		
		phiidCE_MMI(j,i) = phiid_atoms_MMI.str + phiid_atoms_MMI.stx + ...
				phiid_atoms_MMI.sty + phiid_atoms_MMI.sts;
		phiidDC_MMI(j,i) = phiid_atoms_MMI.str + phiid_atoms_MMI.stx + ...
			phiid_atoms_MMI.sty;
		phiidCD_MMI(j,i) = phiidCE_MMI(j,i) - phiidDC_MMI(j,i);
		phiidRed_MMI(j,i) = phiid_atoms_MMI.rtr;
		phiidSyn_MMI(j,i) = phiid_atoms_MMI.sts;
		
		phiid_atoms_CCS = PhiIDFull_Analytical(full_time_lagged_cov, 'CCS');
		
		phiidCE_CCS(j,i) = phiid_atoms_CCS.str + phiid_atoms_CCS.stx + ...
				phiid_atoms_CCS.sty + phiid_atoms_CCS.sts;
		phiidDC_CCS(j,i) = phiid_atoms_CCS.str + phiid_atoms_CCS.stx + ...
			phiid_atoms_CCS.sty;
		phiidCD_CCS(j,i) = phiidCE_CCS(j,i) - phiidDC_CCS(j,i);
		phiidRed_CCS(j,i) = phiid_atoms_CCS.rtr;
		phiidSyn_CCS(j,i) = phiid_atoms_CCS.sts;
		
	end
	% Calculate residuals multi-information (bits)
	J(j,i) = -log2(det(R));
	disp(j)

end 
%}