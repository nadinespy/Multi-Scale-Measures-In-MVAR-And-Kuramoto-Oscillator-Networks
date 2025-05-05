%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bivariate AR(1) model:
%
%    X_t = axx*X_{t-1} + axy*Y_{t-1} + eps_{x,t}
%    Y_t = ayx*X_{t-1} + ayy*Y_{t-1} + eps_{y,t}
%
% Residuals covariance matrix is
%
%         [ sxx  sxy ]
% Sigma = |          |
%         [ sxy  syy ]
%
% RMI is -log(1-kappa^2)/2, where kappa = sxy/sqrt(sxx*syy) is the residuals correlation.
%
% We calculate GC F(Y -> X)
%
% Note that only the axy and ayy AR coefficients enter the expression for F(Y -> X)

defvar('ayy',        0.9    ); % AR coefficient
defvar('axy',       -0.1    ); % AR coefficient
defvar('sxx',        1.0    ); % residuals variance
defvar('syy',        1.0    ); % residuals variance
defvar('rmin',       0.0    ); % rmi range minimum
defvar('rmax',       10.0   ); % rmi range maximum
defvar('rres',       10000  ); % rmi range resolution
defvar('Gnuplot',    true   ); % use Gnuplot? (else Matlab plot)
defvar('pdfviewer', 'mupdf' ); % your favourite PDF viewer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Residuals correlation range

rmi   = linspace(rmin,rmax,rres+1)';
kappa = sqrt(1-exp(-2*rmi));

% F(Y -> X) formula from Supplementary Material for Gutknecht & Barnett, 2023,
% https://academic.oup.com/biomet/article/110/4/933/7035939

sxy = sqrt(sxx*syy)*kappa; % xy component of res. covariance matrix
P   = sxx*(1+ayy^2) - 2*sxy*axy*ayy + syy*(axy^2);
Q   = 2*(sxx*ayy-sxy*axy);
Fxy = log((P+sqrt(P.^2-Q.^2))/2);

% Plot results

figure(1); clf;
plot(rmi, Fxy, 'LineWidth', 2);
xscale log
title({'Granger-Causality (GC) vs. residuals mutual information (RMI)', ...
	'for a bivariate MVAR(1) model'}, 'fontsize', 17, 'interpreter', 'latex');
xlabel('RMI', 'FontSize', 15, 'interpreter', 'latex');
ylabel('$T(Y \to X)$  ', 'rotation', 90, 'FontSize', 15, 'interpreter', 'latex');
set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');  % for tick labels
grid on
	
location = string(strcat(pathout_plots_measures, '2mvar_lag1_lineplot_rmi_gc', '.png'));
exportgraphics(gcf, location);
close all;
