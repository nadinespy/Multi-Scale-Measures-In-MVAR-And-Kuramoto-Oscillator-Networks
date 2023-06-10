%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Need MVGC2

defvar('n',				8    );	% number of variables
defvar('p',				1    );	% time-lags
defvar('w',				1    );	% decay weighting parameter (default is exponential decay of how much past terms are weighted
defvar('rmi',			1    );	% residuals' generalized correlation (multi-information)
defvar('noauto',			false);	% autocorrelation - if 'false', diagonal terms in connectivity matrix are set to zero
defvar('seed',			0    );	% random seed
defvar('coupling_matrices',	false);	% unnormalized connectivity matrices (n x n) for one transition (n x n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed) % seed random number generator

%% calculate integration across range of RMI and C-FACTORS, fixed COUPLING

% choose the n x n coupling matrix with index coup_index from coupling_matrices, 
% and concatenate it as many times as specified by time_lags in third dimension 
% (to account for time-lag)
A = [];
for i = 1:time_lag+1
	A = cat(3,A,coupling_matrices(:,:,coup_index));
end
% specnorm() does take into account decay factor (is it correct that greater time-lag is weighed more?)
A = specnorm(A,1); 

Iv = zeros(length(all_rmi),length(csfac)); % process total variance
Ig = zeros(length(all_rmi),length(csfac)); % process generalised variance (includes covariances)
I = zeros(length(all_rmi),length(csfac));  % integration (process multi-information)
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

	for i = 1:length(csfac)
		C		= var_to_autocov(csfac(i)*A,R,time_lag);	% process covariance matrix (solve DLYAP)
		Iv(j,i)	= sum(log(diag(C(:,:,1))));			% process total variance
		Ig(j,i)	= log(det(C(:,:,1)));				% process generalised variance (includes covariances)
		I(j,i)	= Iv(j,i) - Ig(j,i);				% process multi-information (nats)
		
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
	disp(j);
end
%}
