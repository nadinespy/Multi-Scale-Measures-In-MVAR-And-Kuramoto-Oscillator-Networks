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

% Create a random VAR(1) coefficients matrix

A = randn(n);
A(1:n+1:n^2) = 0;

% Scale A to spectral radius = 1

A = A/max(abs(eig(A)));

% Create a uniform correlation matrix with the same bivariate
% correlation between every pair of variables

R = r*ones(n);
R(1:n+1:n^2) = 1; % set diagonal entries to 1

% Connectivity scale factors in range [0,1) ensures stable VAR(1)

ncsfacs = 1000;
csfac = linspace(0,0.99,ncsfacs)';

% Calculate integration across range of scale factors

I = zeros(ncsfacs,1); % integration (process multi-information)
for i = 1:ncsfacs
	C    = dlyap(csfac(i)*A,R);               % process covariance matrix (solve DLYAP)
	I(i) = sum(log2(diag(C))) - log2(det(C)); % process multi-information (bits) (total correlation)
end

% Calculate residuals multi-information (bits)

J = -log2(det(R));

% Plot integration vs connectivity scaling

figure(1); clf
plot(csfac,I);
yline(J,'r','residuals multi-information');
title(sprintf('Integration vs connectivity scaling (n = %d, r = %g)\n',n,r));
xlabel('connectivity scale factor');
ylabel('integration (bits)');
grid on
