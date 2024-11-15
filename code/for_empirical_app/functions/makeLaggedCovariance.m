function [ lS ] = makeLaggedCovariance(A, c)
  %% Auxiliary function to compute analytically the time-lagged covariance
  % matrix of an AR process with coupling matrix A and noise correlation c
  assert(size(A,1) == size(A,2));
  N = size(A, 1);
  eps = c*ones(N);			% 2*2 matrix (correlations of noise)
  eps(1:(N+1):end) = 1;			% changes diagonal entries to 1.0

  S = dlyap(A, eps);                % solve discrete-time lyapunov equations: ASA' − S + eps = 0 (or S = ASAT + eps); 2*2 matrix; what is S?
  lS = [S, S*A'; A*S, S];           % 4*4 matrix
  lS = 0.5*(lS + lS');              % same 4*4 matrix, should contain the covariances between lagged versions of the time-series (lS considers one time-lag for each variable)
						% overall, I don't know why noise is introduced in this way
						% the rows of lS should be variable 1, variable 1 time lagged, variable 2, variable 2 time lagged, same with columns
                                         
end