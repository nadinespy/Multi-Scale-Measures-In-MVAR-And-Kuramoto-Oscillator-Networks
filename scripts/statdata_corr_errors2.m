% statdata_corr_errors3() obtains time-series data for a Gaussian MVAR(p) process 
% (in a vectorized way) by first getting the time-lagged covariance using connection strength
% and error correlation, extracting the covariance, using it to draw samples of the network,
% and adding the respective A*X_(t-1) term.

function X = statdata_corr_errors2(coupling_matrix, npoints, tau, err);  
	
	settle = 500; %will only keep post-equilibrium data points
	npoints = npoints+settle;
	nvar = size(coupling_matrix,2);
	
	% simulating a time-series using the time-lagged covariance
	S = makeLaggedCovariance(coupling_matrix, err);			    % get time-lagged covariance for a network 
	covariance = S(1:(length(S)/2), 1:(length(S)/2));			    % get same-time covariance of network (by taking the first quadrant of S)
	rng(1);
	X = mvnrnd(zeros([nvar,1])', covariance, npoints)';			 % draw samples from a multi-dimensional correlated Gaussian using
																    % a covariance which already incorporates the connection strengths & error correlation
																    
	% X_t will already entail the error and the covariance of the variables in X 
	% - only the effect of the past term A*X_(t-1) remains to be added.
	
	for t = 1+tau:npoints
		X(:,t) = coupling_matrix*X(:,t-tau) + X(:,t);
	end
	
	X = X(:,settle+1:end);
	
end 