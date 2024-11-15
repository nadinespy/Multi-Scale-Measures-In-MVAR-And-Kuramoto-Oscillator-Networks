function [micro corr_micro] = sim_mvar_model_with_lagged_cov(n_datapoints, full_time_lagged_cov, couplings, ar_model_order)
% sim_mvar_model_with_lagged_cov() simulates a time-series of a multivariate autoregressive 
% network using the full time-lagged covariance

% ...
% -----------------------------------------------------------------------

settle 	= 500;				% keep only post-equilibrium data points
n_datapoints 	= n_datapoints+settle;
nvar 		= size(full_time_lagged_cov,1)/2; 

% *only* the datapoint at t-ar_model_order influences the current datapoint	
covariance 	= full_time_lagged_cov(1:(length(full_time_lagged_cov)/2),1:(length(full_time_lagged_cov)/2));	% get same-time covariance of network 
															% (by taking the first quadrant of S)
rng(1);
micro 		= mvnrnd(zeros([nvar,1])',covariance, n_datapoints)';						% draw samples from a multi-dimensional 
															% correlated Gaussian using a covariance 
															% which already incorporates the connection 
															% strengths & error correlation

% micro_t will already entail the error and the covariance of the variables 
% in micro - only the effect of the past term coupling_matrix*micro_(t-time_lag) remains to be added.

for t = 1+ar_model_order:n_datapoints
	micro(:,t) = couplings*micro(:,t-ar_model_order) + micro(:,t);
end

micro = micro(:, settle+1:end);
% disp(couplings)

corr_micro = corrcoef(micro(1,:),micro(2,:));

end
