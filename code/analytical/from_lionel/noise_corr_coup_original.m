%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Call this script as, e.g.,
%
% >> a = 0.8; b = 0.7; cr = 2; non_monotonic_integration3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Minimal VAR(1) coefficients matrix is A = [a c; 0 b];
%
% Residuals correlation matrix is R = [1 rho; rho 1];

rho = 0:0.01:0.9;

%a = 0.45;
%b = 0.45;
%c = linspace(0.0045,0.45,1001)';

a = 0.8; 
b = 0.7; 
cr = 2;

% Causal connectivity range
c = linspace(-cr,cr,1001)'; 

% Calculate integration for each value of rho across the
% specified range of causal coefficients c

I = minvar1_integration(a,b,c,rho);

% Plot integration vs causal connectivity coefficient

figure(1); clf
plot(c,I);
title(sprintf('Minimal VAR(1): Integration vs causal connectivity (a = %g, b = %g)\n',a,b));
xlabel('causal connectivity coefficient (c)');
ylabel('integration (bits)');
leg = legend(num2str(rho','\\rho = % 3.1f'));
title(leg,'residuals correlation');
grid on

% Function to calculate integration (= process multi-information)

function I = minvar1_integration(a,b,c,rho)

	% Process covariance matrix is C = [p r; r q]
	%
	% Below is the analytic solution of the DLYAP C = A*C*A' + R

	q = 1./(1-b.^2);
	r = (b.*c.*q + rho)./(1-a.*b);
	p = (c.*(2*a.*r +c.*q)+1)./(1-a.^2);

	% Integration is process multi-information (in this case just mutual information)

	I = log2(p)+log2(q)-log2(p.*q-r.^2); % log2 for bits

end
