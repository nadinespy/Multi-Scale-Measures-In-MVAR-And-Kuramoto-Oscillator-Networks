function X = statdata_corr_errors3(A, npoints, err);  
	
	settle = 500; %will only keep post-equilibrium data points
	npoints = npoints+settle;

	nvar = size(A,2);
	S = makeLaggedCovariance(A, err);
	covariance = [S(1,1),S(1,2);S(2,1),S(2,2)];
	X = mvnrnd(zeros([nvar,1])', covariance, npoints)';
	
	X = X(:,settle+1:end);
end 