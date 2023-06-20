function [DD, CE] = get_ce_and_dd(mdim, dim_reduction, cov_matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstrate calculation of Causal Emergence (CE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data (see, e.g., sim_model.m). Do not transform model to
% decorrelated residuals form!
%
% Specify a macroscopic dimension mdim and simulation parameters, or accept
% defaults (see below).
%
% Plotting needs to be commented out, otherwise plotting will happen no matter
% whether gpterm is set to empty or not. (Should fix this at some point.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
defvar('moddir',   tempdir     );  % model directory
defvar('modname', 'sim_model'  );  % model filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('aclmax',   10000       ); % maximum autocovariance lags
defvar('nsamps',   100         ); % number of sample random projections
defvar('iseed',    0           ); % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   []          ); % Gnuplot terminal
defvar('gpscale',  [1,1.2]     ); % Gnuplot scale
defvar('gpfsize',  14          ); % Gnuplot font size
defvar('gpplot',   0           ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

m = mdim;

% Load model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('\n*** loading model from ''%s''... ',modfile);
load(modfile);
fprintf('done\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: all calculations in UNTRANSFORMED coordinates!!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate transfer function

if full_coupling_matrix
	H = var2trfun(full_coupling_matrix, frequency_res);
else
	H = ss2trfun(higher_dim_coupling_matrix, obs_state_var, kalman_var, frequency_res);
end

if full_coupling_matrix
	[G,acl] = var_to_autocov(full_coupling_matrix, noise_corr, aclmax);
else
	[G,acl] = ss_to_autocov(higher_dim_coupling_matrix, obs_state_var, kalman_var, noise_corr, aclmax);
end

fprintf('\n%s: CE calculation (fres = %d, aclags = %d) for m = %d\n',mdescript,frequency_res,acl,m);

fprintf('\nCalculating the Sigma_i ');
st = tic;
[CESRC,CLC] = ac2ces(G);
et = toc(st);
fprintf(' completed in %g seconds\n',et);

% Random projections

rstate = rng_seed(iseed);

if strcmp(dim_reduction, 'grassmanian')
	L = rand_orthonormal(n_nodes,mdim,nsamps); % (orthonormalised) untransformed random linear projections --> take a number of samples
elseif strcmp(dim_reduction, 'pca')			 % do standard PCA
	[coeff, latent, explained] = pcacov(cov_matrix(:,:,1));
end 

rng_restore(rstate);

% Calculate CEs

VRC = chol(noise_corr);
npi = floor(nsamps/10); % for progress indicator

if strcmp(dim_reduction, 'grassmanian')
	CE = zeros(nsamps,1);
	DD = zeros(nsamps,1);
	
	fprintf('\nCalculating causal emergence and dynamical dependence ');
	st = tic;
	for i = 1:nsamps
		[CE(i),DD(i)] = ces2ce(L(:,:,i),H,VRC,CESRC,CLC);
		if ~mod(i,npi), fprintf('.'); end % progress indicator
	end
	et = toc(st);
	fprintf(' completed in %g seconds\n\n',et);

elseif strcmp(dim_reduction, 'pca')
	[CE,DD] = ces2ce(latent,H,VRC,CESRC,CLC);
end
	

% % Plot
%
% gpstem   = fullfile(tempdir,'CE_vs_DD');
% [~,gpname] = fileparts([gpstem '.xxx']); % hack to get fileparts to behave itself
% gp_write(gpstem,[DD CE]);
% gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
% fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
% fprintf(gp,'set title "%s: CE vs DD for m = %d"\n',mdescript,mdim);
% fprintf(gp,'unset key\n');
% fprintf(gp,'set grid\n');
% fprintf(gp,'set xlabel "DD"\n');
% fprintf(gp,'set ylabel "CE" norot\n');
% fprintf(gp,'plot datfile u 1:2 w points pt 7 ps 1 not\n');
% gp_close(gp,gpstem,gpterm,gpplot);


% save plots...?
end
