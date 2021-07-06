function X = statdata_corr_errors2(A, npoints, err);
	
	settle=500; %will only keep post-equilibrium data points
	npoints = npoints+settle;
	
	cov_err = [1 err; err 1];
	nvar = size(A,2);
	tau = floor(size(A,2)/nvar);  
	
	%rng(1);
	%X = normrnd(0,1,nvar,npoints);
	X = zeros([nvar, npoints]);
	
	for t=2:npoints
		X(:,t) = A*X(:,t-tau) + mvnrnd(zeros([1, size(A,1)]), cov_err)';
	end
	
	X = X(:,settle+1:end);
end 

