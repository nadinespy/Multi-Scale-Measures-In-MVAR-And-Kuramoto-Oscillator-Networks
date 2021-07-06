function X = statdata(A, npoints, err)

% -----------------------------------------------------------------------
%   FUNCTION: statdata.m
%   PURPOSE:  Obtain time-series data for a Gaussian MVAR(p) process
%             X_t=A_1*X_{t-1}+A_2*X_{t-2}+...+A_p*X_{t-p}+E_t
%
%   INPUT:  A        -    generalized connectivity matrix. 
%                           A=(A_1 A_2 ... A_p)
%           Omega    -    covariance matrix for E_t
%
%   OUTPUT: X        -    time-series data, rows are variables, columns are
%                         observations
%
%   Adam Barrett May 2010.
% -----------------------------------------------------------------------

settle=500; %will only keep post-equilibrium data points
npoints = npoints+settle;
nvar=size(A,1); 


% simulate time-series of correlated errors (by multiplying M with the upper triangular matrix L obtained from the Cholesky decomposition of the desired correlation matrix R) 
% but: I get different values for corr, which might be due to the fact that chol might fail, if the covarince matrix is singular or near singular. 
%{
mu = 0;
sigma = 1;
M = mu + sigma*randn(N,size(A,2));
R = [1 c; c 1];
L = chol(R);
E = (M*L)';

corr_errors = corrcoef(E(1,:),E(2,:));
%} 

% alternative: simulate errors from correlated multivariate standard normal distribution
% but: this can yield different actual correlations for them; the chosen correlation will only be yielded, if sample size is large 
mu = zeros(1, nvar);
R = [1 err; err 1];

% only necessary for non-standard normal distributions, otherwise R = cov_matrix
standard_dev = ones(1, nvar);								% vector of standard deviations of the errors
cov_matrix = diag(standard_dev)*R*diag(standard_dev);	 % diag() converts std vector to matrix where the diagonal entries are the stds, and all other entries 0
E = mvnrnd(mu,cov_matrix,npoints)';						% simulate correlated errors

corr_errors = corrcoef(E(1,:), E(2,:));						% check correlation

% simulate AR process           
p=floor(size(A,2)/nvar);    

rng(1);
X = normrnd(0,1,nvar,npoints);
%rng(1);
%X_no_err = normrnd(0,1,nvar,npoints);

for i=p+1:npoints
    for j=1:nvar
        for k=1:nvar
            for m=1:p
			%X_no_err(j,i)= X_no_err(j,i)+A(j,k+(m-1)*nvar)*X_no_err(k,i-m);
			X(j,i)= (X(j,i)+A(j,k+(m-1)*nvar)*X(k,i-m))+E(j,i);
            end
        end
    end
end
X = X(:,settle+1:end);

% figure;
% plot(X')
% hold on
% plot(X_no_err')
