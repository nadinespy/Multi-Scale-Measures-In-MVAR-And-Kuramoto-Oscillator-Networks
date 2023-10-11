function [X, Y, psi] = sim_km_oscillators_with_gc_test(coupling_matrix, phase_lag, time_length)
% Default parameters (override on command line - see 'defvar.h')

defvar('C',         tnet5   ); % connectivity matrix, or number of fully-connected oscillators if scalar
defvar('wmean',     0       ); % oscillator frequencies mean
defvar('wsdev',     1/7     ); % oscillator frequencies std. dev.
defvar('wseed',     []      ); % oscillator frequencies random seed (empty for no seeding)
defvar('phase_lag', 0       ); % oscillator phase lag constant
defvar('hseed',     []      ); % oscillator initial phases random seed (empty for no seeding)
defvar('T',         500     ); % simulation time
defvar('Ts',        50      ); % stabilisation time
defvar('dt',        0.01    ); % integration time increment
defvar('nmean',     0.1     ); % oscillator input noise magnitude mean (zero for no noise)
defvar('nsdev',     nmean/3 ); % oscillator input noise magnitude std. dev.
defvar('nseed',     []      ); % oscillator input noise magnitude random seed (empty for no seeding)
defvar('Iseed',     []      ); % oscillator input noise random seed (empty for no seeding)
defvar('smode',    'Euler'  ); % simulation mode: 'Euler' or 'RK4'
defvar('hifreq',    0.01    ); % hi-pass frequency, or 0 for no high-pass
defvar('pdto',      1       ); % polynomial detrend order
defvar('dsf',       2       ); % downsample factor
defvar('varmomax',  40      ); % maximum VAR model order for model order selection
defvar('varmosel', 'HQC'    ); % model order selection critierion: 'AIC', 'BIC', 'HQC', or 'LRT'
defvar('stest',    'F'      ); % statistical test: 'LR' or 'F'
defvar('mhtc',     'FDRD'   ); % multiple-hypothesis test correction
defvar('sig_level', 0.05    ); % significance level
defvar('fignum',    1       ); % Matlab figure number

% Connectivity

C = length(coupling_matrix);
T = time_length;
if isscalar(C)  % fully connected
	N = C;
	C = ones(N);
else
	N = size(C,1);
	assert(size(C,2) == N,'Connectivity matrix must be square');
end
C(1:N+1:N^2) = 0; % zeros on diagonal (no self-connections!)

% Random Kuramoto parameters

if ~isempty(wseed), rstate = rng(wseed); end
w = wmean + wsdev*randn(N,1); % oscillator frequencies normally distributed with mean wmean and std. dev wsdev
if ~isempty(wseed), rng(rstate); end

if ~isempty(hseed), rstate = rng(hseed); end
h0 = pi*(2*rand(N,1)-1);      % initial phases uniform on [-pi,pi]
if ~isempty(hseed), rng(rstate); end

ns = round(Ts/dt);
assert(ns > 0,'Stabilisation time too short, or time increment too large!');
na = round((T+Ts)/dt);
n  = na-ns;
T  = n*dt;  % adjusted simulation time
Ts = ns*dt; % adjusted total time
Ta = na*dt; % adjusted stabilisation time

if nmean > 0 % with input noise
	if ~isempty(nseed), rstate = rng(nseed); end
	lnv = log(1+nsdev^2/nmean^2);
	nmag = (lognrnd(log(nmean)-lnv/2,sqrt(lnv),1,N));
	if ~isempty(nseed), rng(rstate); end
	if ~isempty(Iseed), rstate = rng(Iseed); end
	I = nmag'.*randn(N,n);
	if ~isempty(Iseed), rng(rstate); end
else
	I = []; % no input
end

% Run Kuramoto simulation with specified parameters

fprintf('\n');

st = tic;
%[h,r] = kuramoto(N,na,dt,w,K/N,a,h0,I,smode);
[h,r, psi] = kuramoto(N,w,coupling_matrix,phase_lag,n,dt,I, 'Euler');

et = toc(st);
fprintf('%s method : %g seconds\n',smode,et);

% Truncate and transpose

t = linspace(0,T,n-ns)';
r = r(ns+1:end)';
h = h(:,ns+1:end)';
psi = psi(:,ns+1:end)';

% Polynomial detrend and normalise by variance

fprintf('\nPolynomial detrend at order %d\n',pdto);
y = detrend(h,pdto);
y = demean(y',true)'; %
pdstr = sprintf('Order %d polynomial detrend',pdto);

% Optionally high-pass

if hifreq > 0
	fprintf('\nHigh-pass filter at %g/2pi\n',hifreq);
	hfstr = sprintf('High-pass filter at %g/2\\pi',hifreq);
	[fb,fa] = butter(2,hifreq/pi,'high');
	x = filtfilt(fb,fa,y);
else
	x = y;
end

% Plots

figure(fignum); clf
sgtitle(sprintf('%s simulation\n',smode));

subplot(3,1,1)
plot(t,r);
title('Order parameter magnitude');
xlabel('time');
ylabel('r');
ylim([0,1]);

subplot(3,1,2)
plot(t,y);
title(pdstr);
xlabel('time');
ylabel('phase');
ymax = 1.05*max(abs(y(:)));
ylim([-ymax,ymax]);

if hifreq > 0
	subplot(3,1,3)
	plot(t,x);
	title(hfstr);
	xlabel('time');
	ylabel('phase');
	xmax = 1.05*max(abs(x(:)));
	ylim([-xmax,xmax]);
end

% Downsample and normalise by variance (again)

fprintf('\nDownsample by a factor of %d\n',dsf);
m = floor((n-ns)/dsf);
X = zeros(N,m);
for i = 1:N
	X(i,:) = downsample(x(:,i),dsf);
end
X = demean(X,true);

%Y = zeros(size(r,2),m);
Y = downsample(r,dsf);

% VAR model order selection

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,varmomax,'LWR',[],[],2);
switch varmosel
	case 'AIC', varmo = moaic;
	case 'BIC', varmo = mobic;
	case 'HQC', varmo = mohqc;
	case 'LRT', varmo = molrt;
	otherwise, error('Unknown VAR model order selection criterion');
end
fprintf('\nUsing VAR model order (%s) = %d\n',varmosel,varmo);

% VAR modelling

[A,V] = tsdata_to_var(X,varmo,'OLS');
rho = specnorm(A);
fprintf('\nSpectral norm = %g\n',rho);

% Granger causality calculation (time domain)

F = var_to_pwcgc(A,V);
fprintf('\nPairwise-conditional GC =\n\n'); disp(F)

% Statistical inference

tstat = var_to_pwcgc_tstat(X,V,varmo,'OLS',stest);
pval  = mvgc_pval(tstat,stest,1,1,N-2,varmo,m,1); % for pairwise-conditional, nx = 1, ny = 1, nz = nvars-2
sig   = significance(pval,sig_level,mhtc);

pdata = {abs(coupling_matrix),F,sig};
ptitle = {'Connectivity','PWCGC','Significant'};
maxp = [nanmax(coupling_matrix(:)),nanmax(F(:)),1];
plot_gc(pdata,ptitle,[],maxp,3,[0.6,2.5]);

end
