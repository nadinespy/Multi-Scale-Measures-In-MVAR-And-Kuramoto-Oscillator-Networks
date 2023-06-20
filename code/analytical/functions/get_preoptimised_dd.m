function get_preoptimised_dd(mdim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence pre-optimisation (via proxy DD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data (see, e.g., sim_model.m).
%
% Specify a macroscopic dimension mdim and optimisation parameters, or accept
% defaults (see below). After running this script, run optimise_dd.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',    tempdir    );  % model directory
defvar('modname',  'sim_model' );  % model filename root
defvar('poptdir',   tempdir    );  % pre-optimisation directory
defvar('poptname', 'preopt_dd' );  % pre-optimisation filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('iseed',     0          ); % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('nrunsp',    100        ); % pre-optimisation runs (restarts)
defvar('nitersp',   10000      ); % pre-optimisation iterations
defvar('gdesp',     2          ); % gradient-descent ES version (1 or 2)
defvar('gdsig0p',   1          ); % pre-optimisation (gradient descent) initial step size
defvar('gdlsp',     2          ); % gradient-descent "line search" parameters
defvar('gdtolp',    1e-10      ); % gradient descent convergence tolerance
defvar('histp',     true       ); % calculate optimisation history?
defvar('ppp',       false      ); % parallelise multiple runs?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',    []         ); % Gnuplot terminal (empty for no plots. 'x11' for plotting)
defvar('gpscale',   [Inf,0.6]  ); % Gnuplot scale
defvar('gpfsize',   14         ); % Gnuplot font size
defvar('gpplot',    0          ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

m = mdim;

% Load model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('\n*** loading model from ''%s''... ',modfile);
load(modfile);
fprintf('done\n\n');

% Transform model to decorrelated and normalised form, and calculate CAK sequence

if full_coupling_matrix
	[decorr_coupling_matrix, decorr_noise] = transform_var(full_coupling_matrix,noise_corr);       % transform model to decorrelated-residuals form
	[decorr_higher_dim_coupling_matrix, decorr_obs_state_var, decorr_kalman_var] = var_to_ss(decorr_coupling_matrix);               % equivalent ISS model
	CAK = decorr_coupling_matrix;
else
	[decorr_higher_dim_coupling_matrix, ...
		decorr_obs_state_var, ...
		decorr_kalman_var, ...
		decorr_noise] = transform_ss(higher_dim_coupling_matrix, ...
		observed_state_var, ...
		kalman_var, ...
		noise_corr);  % transform model to decorrelated-residuals form
	
	CAK = iss2cak(decorr_higher_dim_coupling_matrix, decorr_obs_state_var, decorr_kalman_var);
end

fprintf('%s: pre-optimisation for m = %d\n\n',mdescript,m);

% Initialise optimisations

rstate = rng_seed(iseed);
L0p = rand_orthonormal(n_nodes,m,nrunsp); % initial (orthonormalised) random linear projections
rng_restore(rstate);

% Multiple optimisation runs

st = tic;
[doptp,Lp,convp,ioptp,soptp,cputp,ohistp] = opt_gd_ddx_mruns(CAK,L0p,nitersp,gdesp,gdsig0p,gdlsp,gdtolp,histp,ppp);
et = toc(st);

% Inverse-transform Lp back for un-decorrelated residuals

Loptp = itransform_subspace(Lp,noise_corr);

% print runs
fprintf('\noptimal dynamical dependence =\n'); disp(doptp');
fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));
fprintf('CPU secs per run = %7.4f +- %6.4f\n\n',mean(cputp),std(cputp));

% Plot optimisation histories

if histp && ~isempty(gpterm)
	gptitle  = sprintf('Pre-optimisation history: %s, m = %d',mdescript,m);
	gpstem   = fullfile(tempdir,'preopt_hist');
	gp_opthist({ohistp},nitersp,true,true,{'Pre-optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

goptp = gmetrics(Loptp);
if ~isempty(gpterm)
	gptitle = sprintf('Inter-preoptimum distance: %s, m = %d',mdescript,m);
	gpstem = fullfile(tempdir,'preopt_iodist');
	gp_iodist(goptp,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);
end

% Save pre-optimisation results

clear n m st et tcpu rstate gpplot gpterm gpscale gpfsize gptitle gpstem
if histp
	poptfile = fullfile(poptdir,[poptname '_mdim_' num2str(mdim) '_H.mat']);
else
	poptfile = fullfile(poptdir,[poptname '_mdim_' num2str(mdim) '_N.mat']);
end

fprintf('\n*** saving pre-optimisation results in ''%s''... ',poptfile);
save(poptfile);
fprintf('done\n\n');

end
