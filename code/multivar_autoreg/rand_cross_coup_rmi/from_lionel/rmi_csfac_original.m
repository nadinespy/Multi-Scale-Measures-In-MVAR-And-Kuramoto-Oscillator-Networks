%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Need MVGC2

defvar('n',      10   );
defvar('p',      4    );
defvar('w',      1    );
defvar('rmi',    1    );
defvar('noauto', false);
defvar('seed',   0    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed) % seed random number generator

% Create a random VAR(1) coefficients matrix with spectral radius = 1

A = var_rand(n,p,1,w);
if noauto
	d = 1:n+1:n^2;
	for k = 1:p
		Ak = A(:,:,k);
		Ak(d) = 0;
		A(:,:,k) = Ak;
	end
end
A = specnorm(A,1);

% Create a uniform correlation matrix with the same bivariate
% correlation between every pair of variables

R = corr_rand(n,rmi);

% Connectivity scale factors in range [0,1) ensures stable VAR(1)

ncsfacs = 1000;
csfac = linspace(0,0.99,ncsfacs)';

% Calculate integration across range of scale factors

Iv = zeros(ncsfacs,1); % process total variance
Ig = zeros(ncsfacs,1); % process generalised variance (includes covariances)
I = zeros(ncsfacs,1);  % integration (process multi-information)
for i = 1:ncsfacs
	C    = var_to_autocov(csfac(i)*A,R,0);  % process covariance matrix (solve DLYAP)
	Iv(i) = sum(log(diag(C)));              % process total variance
	Ig(i) = log(det(C));                    % process generalised variance (includes covariances)
	I(i)  = Iv(i) - Ig(i);                  % process multi-information (nats)
end

% Plot integration, etc. vs connectivity scaling

figure(1); clf
plot(csfac,[I Iv Ig]);
yline(rmi,'g','residuals multi-information');
title(sprintf('Integration vs connectivity scaling (n = %d, model order = %d, rmi = %g)\n',n,p,rmi));
legend({'integration','total variance','generalised variance'});
xlabel('connectivity scale factor');
ylabel('integration (nats)');

grid on
