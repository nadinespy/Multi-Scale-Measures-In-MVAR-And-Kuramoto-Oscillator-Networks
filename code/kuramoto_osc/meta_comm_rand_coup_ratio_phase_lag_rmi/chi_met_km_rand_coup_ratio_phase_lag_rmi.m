%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mec_km_rand_cou_ratio_phase_lag_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this script, run those scripts before (in that order):
% >> params_dd_ce_co_info_rand_phase_lag_rmi.m;
% >> mec_metacommkm_rand_phase_lag_rmi_local_dirs.m
%
% ADD DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DISPLAY PARAMETERS

fprintf('\nParameters:\n-----------\n\n');
fprintf('current directory: %s\n',		pwd);
fprintf('base directory: %s\n',			basedir);
fprintf('results_filename: %s\n',		results_filename);
fprintf('measure: %s\n',			measure);
fprintf('ar_model_order: %d\n',			ar_model_order);
fprintf('seed: %d\n',				seed);
% -------------------------------------------------------------------------
fprintf('n_nodes: %d\n',			n_nodes);
fprintf('n_comms: %d\n',			n_comms);
fprintf('n_inter_comm_coups: %d\n',		n_inter_comm_coups);
fprintf('n_nodes_per_comm: %d\n',		n_nodes_per_comm);
% -------------------------------------------------------------------------
fprintf('phase_lag_range: %g - %g\n',		phase_lag_range(1), ...
	phase_lag_range(2));
fprintf('n_phase_lag: %d\n',			n_phase_lag);
% -------------------------------------------------------------------------
fprintf('rmi_range: %g - %g\n',			rmi_range(1), rmi_range(2));
fprintf('n_rmi: %d\n',				n_rmi);
fprintf('n_samples_noise_corr: %d\n',		n_samples_noise_corr);
fprintf('noise_mag: %d\n',			noise_mag);
fprintf('adjust_rmi: %d\n',			adjust_rmi);
% -------------------------------------------------------------------------
fprintf('coup_ratio_mean_range: %g - %g\n',	coup_ratio_mean_range(1), ...
	coup_ratio_mean_range(2));
fprintf('n_coup_ratio_mean: %d\n',		n_coup_ratio_mean);
fprintf('n_samples_coup: %d\n',			n_samples_coup);
fprintf('coup_ratio_mean_rel_dev: %d\n',	coup_ratio_mean_rel_dev);
fprintf('global_coup: %d\n',			global_coup);
% -------------------------------------------------------------------------
fprintf('n_sim_points: %d\n',			n_sim_points);
fprintf('n_equil_points: %d\n',			n_equil_points);
fprintf('dt_sim_points: %d\n',			dt_sim_points);
fprintf('sim_mode: %d\n',			sim_mode);
% -------------------------------------------------------------------------
fprintf('n_local_gc: %d\n',			n_local_gc);
fprintf('emp_sample: %d\n',			emp_sample);
fprintf('stabilise_care: %d\n',			stabilise_care);
% -------------------------------------------------------------------------
fprintf('osc_freq_mean: %d\n',			osc_freq_mean);
fprintf('osc_freq_rel_dev: %d\n',		osc_freq_rel_dev);
% -------------------------------------------------------------------------
fprintf('\n');

% --------------------------------------------
% measure-specific parameters
% --------------------------------------------

%% SET PATHS

% run kvar (also calls SSDI, MVGC2, GVMAT & GPMAT startup)
run(fullfile(basedir, 'packages_and_code_repos', 'kvar', 'startup'));

addpath(fullfile(basedir, 'functions'));
addpath(fullfile(basedir, 'scripts', 'kuramoto_osc', ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi'));

% {
javaaddpath(fullfile(basedir, 'scripts', 'kuramoto_osc', ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi', 'infodynamics_with_iit.jar'));
javaaddpath(fullfile(basedir, 'scripts', 'kuramoto_osc', ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi', 'commons-math3-3.5.jar'));
javaaddpath(fullfile(basedir, 'functions', 'infodynamics_with_iit.jar'));
javaaddpath(fullfile(basedir, 'functions', 'commons-math3-3.5.jar'));

try
	javaaddpath(fullfile(basedir, 'scripts', 'kuramoto_osc', ...
	'meta_comm_rand_coup_ratio_phase_lag_rmi', 'infodynamics_with_iit.jar'));
	disp('Path added successfully.');
catch exception
	disp(['Error adding path: ' exception.message]);
end

javaClassPath = javaclasspath;
disp(javaClassPath);
%}

rng(seed); % seed random number generator

%% MEMORY PRE-ALLOCATION

% coupling/noise zero matrices
coup_noise_corr_zeros 		= zeros(n_samples_coup, n_samples_noise_corr);
blank_struct.global_coup	= coup_noise_corr_zeros;
blank_struct.chimera_index	= coup_noise_corr_zeros;
blank_struct.metastability	= coup_noise_corr_zeros;

results = repmat(blank_struct, n_coup_ratio_mean, n_phase_lag, n_rmi); % preallocate

w = whos('results');
fprintf('Pre-allocated size of ''results'' = %.4f MB\n', ...
	w.bytes/1000/1000);

%% SET UP KM OSCILLATOR SYSTEM

seed_osc_freq		= 0;	% oscillator frequencies random seed (zero for no seeding)
seed_coup_ratio_mean	= 0;	% Shanahan connectivity parameter random seed (zero for no seeding)
seed_noise_corr	= 0;	% oscillator input noise correlation random seed (zero for no seeding)
seed_noise_gen		= 0;	% oscillator input noise generation random seed (zero for no seeding)
seed_init_phase	= 0;	% oscillator initial phases random seed (zero for no seeding)

all_rmi = 2.^linspace(rmi_range(1), rmi_range(2), n_rmi);		% rmi values (exponential scale)
all_phase_lag = flip(pi/2 - linspace(phase_lag_range(1), ...		% phase lag, flip order so that values increase
	phase_lag_range(2), n_phase_lag));
all_coup_ratio_mean = linspace(coup_ratio_mean_range(1), ...		% coupling ratio mean values
	coup_ratio_mean_range(2), n_coup_ratio_mean);

osc_freq_dev = osc_freq_rel_dev * osc_freq_mean;			% oscillator frequencies deviation

% oscillator natural frequencies
if osc_freq_dev > 0
	seed_state = rng_seed(seed_osc_freq);
	osc_freq = pi*ms_betarnd(osc_freq_mean/pi, ...			% oscillator frequencies Beta
		osc_freq_dev, n_nodes, 1);				% distributed on [0,pi] with
	rng_restore(seed_state);					% mean osc_freq_mean and
else									% std. dev osc_freq_dev
	osc_freq = osc_freq_mean * ones(n_nodes, 1);
end
sampling_freq  = 1/dt_sim_points;					% sampling frequency
nyquist_freq = sampling_freq/2;						% Nyqvist frequency
freq_hz = nyquist_freq * osc_freq/pi;					% frequencies in Hz

% simulation times
dt_equil_points = floor(n_equil_points / dt_sim_points);		% equilibriation time increments
assert(dt_equil_points >= 0, ...
	'Equilibriation time too short, or time increment too large!');
n_equil_points = dt_equil_points * dt_sim_points;			% adjusted equilibriation time

sim_time_increments = floor(n_sim_points/dt_sim_points);		% simulation time increments
assert(sim_time_increments > 0, ...
	'simulation time too short, or time increment too large!');
n_sim_points = sim_time_increments * dt_sim_points;			% adjusted simulation time

tot_sim_time_increments = dt_equil_points + sim_time_increments;	% total simulation time increments

%% MODEL PARAMETER SWEEP (RMI, COUPLING RATIO MEAN, & PHASE LAG)

for h = 1:n_coup_ratio_mean;

	% set up community-structured weighted network as in Shanahan (2009),
	% if coupling ratio has zero variance (i.e., all inter-community
	% couplings will be the same, and all intra-community couplings
	% will be the same); this implies that only one particular coupling
	% matrix will be used, with a number of different noise correlation
	% matrices

	% coupling ratio mean deviation
	coup_ratio_mean_dev = coup_ratio_mean_rel_dev * all_coup_ratio_mean(h);

	if coup_ratio_mean_rel_dev == 0 && n_samples_coup == 1
		seed_state = rng_seed(seed_coup_ratio_mean);
		[coupling_matrix, communities] = shanahan_network(n_comms, ...
			n_nodes_per_comm, ...
			n_inter_comm_coups, ...
			all_coup_ratio_mean(h), ...
			coup_ratio_mean_dev, ...
			global_coup);
		rng_restore(seed_state);
	end

	for i = 1:length(all_phase_lag);

		for j = 1:n_rmi;
			fprintf('\nCOUPLING RATIO MEAN %d of %d = %g\n', h, n_coup_ratio_mean, ...
				all_coup_ratio_mean(h));
			fprintf('\nPHASE LAG %d of %d = %g\n', i, n_phase_lag, ...
				all_phase_lag(i));
			fprintf('\n\tRMI %d of %d = %g\n', j, n_rmi, all_rmi(j));
			st = tic;

			results(h,i,j).rmi_for_this_sim		= all_rmi(j);
			results(h,i,j).coup_ratio_mean_for_this_sim	= all_coup_ratio_mean(h);
			results(h,i,j).phase_lag_for_this_sim		= all_phase_lag(i);

			for k = 1:n_samples_coup;
				fprintf('\t\tcouplings sample %d of %d\n', ...
					k, n_samples_coup);

				% set up community-structured weighted network as in Shanahan (2009),
				% if coupling ratio has non-zero variance (i.e., inter-community
				% couplings will have variance around some mean, and all
				% intra-community couplings will have variance around some mean);
				% this implies we'll generate multiple different coupling matrices,
				% with the same coupling ratio mean and variance, as well as a
				% number of different noise correlation matrices
				if coup_ratio_mean_rel_dev ~= 0 && n_samples_coup > 1
					seed_state = rng_seed(seed_coup_ratio_mean);
					[coupling_matrix, communities] = shanahan_network(n_comms, ...
						n_nodes_per_comm, ...
						n_inter_comm_coups, ...
						all_coup_ratio_mean(h), ...
						coup_ratio_mean_dev, ...
						global_coup);
					rng_restore(seed_state);
				end
				% make sure coupling matrix has been generated
				error_msg = ['Coupling matrix does not exist.', ...
					' "coup_ratio_mean_rel_dev" might be 0, while "n_samples_coup" > 1,', ...
					' or "coup_ratio_mean_rel_dev" is > 0, while "n_samples_coup" is 1.'];

				assert(exist('coupling_matrix', 'var'), error_msg);

				for p = 1:n_samples_noise_corr;
					fprintf('\t\t\tnoise corr sample %d of %d\n', ...
						p, n_samples_noise_corr);

					%------------------------------------------------------------------------
					% COUPLING & NOISE CORRELATION MATRICES
					%------------------------------------------------------------------------

					% generate noise covariance matrix
					seed_state = rng_seed(seed_noise_corr);
					if adjust_rmi
						noise_cov = sqrt(noise_mag) * corr_rand(n_nodes, -all_rmi(j));
					else
						noise_cov = sqrt(noise_mag) * corr_rand(n_nodes, +all_rmi(j));
					end
					rng_restore(seed_state);

					% generate correlated Gaussian white noise
					seed_state = rng_seed(seed_noise_gen);
					corr_noise_time_series = randn(tot_sim_time_increments, n_nodes) * chol(noise_cov);
					rng_restore(seed_state);

					% initial phases uniform on [-pi,pi]
					seed_state = rng_seed(seed_init_phase);
					initial_phases = pi*(2*rand(1, n_nodes)-1);
					rng_restore(seed_state);

					% run Kuramoto simulations with specified parameters
					[phases, order_param] = kuramoto(n_nodes, ...
						tot_sim_time_increments, ...
						dt_sim_points, osc_freq, ...
						coupling_matrix, ...
						all_phase_lag(i), ...
						initial_phases, ...
						corr_noise_time_series, ...
						sim_mode);

					% truncate equilibriation
					phases = phases(:, dt_equil_points+1:end);
					order_param = order_param(dt_equil_points+1:end);

					% community order parameters (synchronisation)
					% --> macroscopic variable
					order_param_comms = zeros(n_comms, sim_time_increments);
					for o = 1:n_comms
						order_param_comms(o,:) = hypot(mean(cos(phases(...
							communities{o},:))), mean(sin(phases(communities{o},:))));
					end
					
					% calculate global metastability:
					% mean of variance over time for each synchrony
					results(h,i,j).global_metastability(k,p) = mean(var(order_param_comms')); % community_sync': time x synchronies
				
					% calculate global chimera index
					% mean of variance over synchronies for each time-point
					results(h,i,j).global_chimera_index(k,p) = mean(var(order_param_comms));  % community_sync: synchronies x time

					% global coupling
					results(h,i,j).global_coup(k,p) = sqrt(mean(...
						coupling_matrix(diag(true(n_nodes, 1))).^2)/ ...
						mean(coupling_matrix(:).^2));
					
					fprintf('\n\tglobal_chimera_index = %g\n', results(h,i,j).global_chimera_index(k,p));
					fprintf('\n\tglobal_metastability = %g\n', results(h,i,j).global_metastability(k,p));
					fprintf('\n\tglobal_coup = %g\n', results(h,i,j).global_coup(k,p));

					clear noise_cov;
					clear phases;
					clear coupling_matrix;
				end
			end

			w = whos('results');
			et = toc(st);
			fprintf('\t\tcoupling/noise loop completed in %g seconds; size of ''results'' = %.4f MB\n', ...
				et, w.bytes/1000/1000);
		end
	end
end

% save results
pathout_data_measures  = fullfile(basedir, 'results', 'analyses', ...
        'kuramoto_osc', 'meta_comm_rand_coup_ratio_phase_lag_rmi', model_specific_path);

% save results
results_filepath = fullfile(pathout_data_measures, ...
	[model_specific_filename '_' results_filename '.mat']);
fprintf('\nSaving results file ''%s'' ...', ...
	results_filepath);
save(results_filepath, 'results', '-v7.3');
fprintf(' one\n');



