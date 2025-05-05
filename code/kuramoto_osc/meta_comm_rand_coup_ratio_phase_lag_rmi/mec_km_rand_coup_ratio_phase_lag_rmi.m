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
coup_noise_corr_zeros = zeros(n_samples_coup, n_samples_noise_corr);
noise_corr_zeros = zeros(n_nodes, n_nodes, n_samples_coup, ...
        n_samples_noise_corr);
coup_zeros = zeros(n_nodes, n_nodes, ar_model_order, ...
        n_samples_coup, n_samples_noise_corr);

% the same for all measures
blank_struct.model				= model_specific_filename;
blank_struct.n_nodes				= n_nodes;
blank_struct.n_comms				= n_comms;
blank_struct.n_inter_comm_coups			= n_inter_comm_coups;
blank_struct.n_nodes_per_comm			= n_nodes_per_comm;
blank_struct.ar_model_order			= ar_model_order;
blank_struct.rmi_for_this_sim			= 0;
blank_struct.coup_ratio_mean_for_this_sim	= 0;
blank_struct.phase_lag_for_this_sim		= 0;

blank_struct.phase_lag_range			= phase_lag_range;

blank_struct.rmi_range				= rmi_range;
blank_struct.noise_mag				= noise_mag;
blank_struct.noise_corr				= noise_corr_zeros;

blank_struct.coup_ratio_mean_range		= coup_ratio_mean_range;
blank_struct.coup_ratio_mean_rel_dev		= coup_ratio_mean_rel_dev;
blank_struct.global_coup			= global_coup;
blank_struct.coupling_matrices			= coup_zeros;

blank_struct.n_sim_points			= n_sim_points;
blank_struct.n_equil_points			= n_equil_points;
blank_struct.dt_sim_points			= dt_sim_points;
blank_struct.sim_mode				= sim_mode;

blank_struct.n_local_gc				= n_local_gc;
blank_struct.emp_sample				= emp_sample;
blank_struct.stabilise_care			= stabilise_care;

blank_struct.osc_freq_mean			= osc_freq_mean;
blank_struct.osc_freq_rel_dev			= osc_freq_rel_dev;

blank_struct.seed				= seed;

% measure-dependent
switch measure

	case 'multi_info'
		blank_struct.MultiInfo = coup_noise_corr_zeros;

	case 'dd_ce_co_info'
		if any(strcmp(dim_reduction,'comm_sync'))
			blank_struct.m_dim_ce_dd_co_info	= n_comms;
			blank_struct.DDCommSync			= coup_noise_corr_zeros;
			blank_struct.ShannonCECommSync		= coup_noise_corr_zeros;
			blank_struct.CoInfoCommSync		= coup_noise_corr_zeros;
		else
			error('Undefined dim_reduction, must be "comm_sync".')
		end

	case 'phiid_measures_mmi'
		blank_struct.poly_detrend_order		= poly_detrend_order;
		blank_struct.high_pass_freq		= high_pass_freq;
		blank_struct.n_features			= n_features;
		blank_struct.PhiID_DoubleRed_MMI	= coup_noise_corr_zeros;
		blank_struct.PhiID_CE_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_DC_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_CD_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_UC_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_Syn_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_Transfer_MMI		= coup_noise_corr_zeros;

	case 'phiid_measures_ccs'
		blank_struct.poly_detrend_order		= poly_detrend_order;
		blank_struct.high_pass_freq		= high_pass_freq;
		blank_struct.PhiID_DoubleRed_CCS	= coup_noise_corr_zeros;
		blank_struct.PhiID_CE_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_DC_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_CD_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_UC_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_Syn_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_Transfer_CCS	= coup_noise_corr_zeros;

	case 'integrated_info_measures'
		if any(strcmp(list_integrated_info_measures, ...
				'AverageCorrelation'   )), ...
				blank_struct.AverageCorr		= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'IntegratedInformation')), ...
				blank_struct.IntegratedInfo		= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'IntegratedInteraction')), ...
				blank_struct.IntegratedInteraction	= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'DecoderIntegration'   )), ...
				blank_struct.DecoderIntegration		= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'CausalDensity'        )), ...
				blank_struct.CausalDensity		= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'IntegratedSynergy'    )), ...
				blank_struct.IntegratedSynergy		= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'TimeDelayedMutualInfo')), ...
				blank_struct.TimeDelayedMI		= coup_noise_corr_zeros; end

	otherwise error('Unknown measure ''%s''', measure);
end % switch measure

results = repmat(blank_struct, n_coup_ratio_mean, n_phase_lag, n_rmi); % preallocate

w = whos('results');
fprintf('Pre-allocated size of ''results'' = %.4f MB\n', ...
	w.bytes/1000/1000);

%% SET UP KM OSCILLATOR SYSTEM

seed_osc_freq		= 0;	% oscillator frequencies random seed (zero for no seeding)
seed_coup_ratio_mean	= 0;	% Shanahan connectivity parameter random seed (zero for no seeding)
seed_noise_corr		= 0;	% oscillator input noise correlation random seed (zero for no seeding)
seed_noise_gen		= 0;	% oscillator input noise generation random seed (zero for no seeding)
seed_init_phase		= 0;	% oscillator initial phases random seed (zero for no seeding)

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

			results(h,i,j).rmi_for_this_sim			= all_rmi(j);
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

					% compute time-lagged covariance matrix of KM time series
					% phases_autocov_sequence = tsdata_to_autocov(phases,1);

					results(h,i,j).noise_corr(:,:,k,p) = noise_cov;
					results(h,i,j).coupling_matrices(:,:,:,k,p) = coupling_matrix;

					% HOW TO GET GLOBAL COUPLING?
					%results(h,i,j).global_coup(k,p) = sqrt(mean(...
					%	coupling_matrix(diag(true(n_nodes, 1))).^2)/ ...
					%	mean(coupling_matrix(:).^2));

				switch measure

					%------------------------------------------------------------------------
					% MULTI-INFORMATION
					%------------------------------------------------------------------------
					case 'multi_info'

						% process multi-information (nats)
						results(h,i,j).MultiInfo(k,p) = sum(log(diag(...
							phases_autocov_sequence(:,:,1)))) - ...
							logdet(phases_autocov_sequence(:,:,1));

						fprintf('\n\tMultiInfo = %g\n', results(h,i,j).MultiInfo(k,p));

					%------------------------------------------------------------------------
					% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
					%------------------------------------------------------------------------
					case 'dd_ce_co_info'

						for z = 1:length(dim_reduction)

							if strcmp(dim_reduction{z}, 'comm_sync')

								DD = nan(n_local_gc, 1);
								sample_counter = 0;
								n_fail = 0;
								while sample_counter < n_local_gc
									if emp_sample
										phases_sample = phases(:, ...
											randi(sim_time_increments));
									else
										phases_sample = pi*(2 * rand(n_nodes,1)-1);
									end

									% locally linearised Kuramoto process
									local_lin_phases = kosgrad(phases_sample, ...
										coupling_matrix, all_phase_lag(i));
									local_stability_index = max(real(eig(...
										local_lin_phases)));

									% locally linearised macro variable
									local_lin_order_param_comms = kosmacrograd(...
										phases_sample, ...
										communities);

									% compute Dynamical Dependence (DD) for locally
									% linearised micro & macro variables
									[DD_temp, err] = vou_to_dd(local_lin_phases, ...
										noise_cov, ...
										local_lin_order_param_comms, ...
										stabilise_care);

									if err == 0 % success!
										sample_counter = sample_counter+1;
										DD(sample_counter) = DD_temp;
										if sample_counter == n_local_gc, break;
										end
									else
										n_fail = n_fail+1;
									end
								end

								DD_mean = mean(DD);
								DD_std = std(DD);

								results(h,i,j).DDCommSync(k,p)= DD_mean;
								results(h,i,j).local_stability_index(:,:,:,k,p) = local_stability_index;

								fprintf('\n\tDDCommSync = %g\n', ...
									results(h,i,j).DDCommSync(k,p));
							end
						end

					%------------------------------------------------------------------------
                              % PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, CAUSAL DECOUPLING
                              %------------------------------------------------------------------------
					case {'phiid_measures_mmi', 'phiid_measures_ccs', ...
							'integrated_info_measures'}

						% do polynomial detrending and normalise by variance to make KM
						% time-series Gaussian

						phases_detrended = detrend(phases, poly_detrend_order);
						phases_detrended = demean(phases_detrended',true)';

						% Optionally high-pass

						if high_pass_freq > 0
							[fb,fa] = butter(2, high_pass_freq/pi,'high');
							phases_transformed = filtfilt(fb, fa, phases_detrended);
						else
							phases_transformed = phases_detrended;
						end

						if strcmp(measure, 'phiid_measures_mmi') || ...
								strcmp(measure,'phiid_measures_ccs')

							% reduce dimensionality of system to 8, as PhiID-based measures
							% can't be computed for large systems
							n_features = 8;
							reconstruction_ica = rica(phases_transformed', n_features, ...
								'IterationLimit', 1000);
							phases_rica = (phases_transformed' * ...
								reconstruction_ica.TransformWeights)';

							% compute time-lagged covariance matrix of transformed KM time series
							phases_autocov_sequence = tsdata_to_autocov(phases_rica,1);

							% construct full time-lagged covariance matrix (block-toeplitz form)
							phases_lag01_autocov = [phases_autocov_sequence(:,:,1), ...
								phases_autocov_sequence(:,:,2)'; ...
								phases_autocov_sequence(:,:,2), ...
								phases_autocov_sequence(:,:,1)];

							if contains(measure, 'mmi')
								red_func = 'MMI';
							elseif contains(measure, 'ccs')
								red_func = 'CCS';
							else
								error('No redundancy function defined.');
							end

							PhiID_Atoms		= PhiIDFullAnalytical(phases_lag01_autocov, red_func);

							PhiID_CE		= PhiID_Atoms.str + PhiID_Atoms.stx + ...
								PhiID_Atoms.sty + PhiID_Atoms.sts;
							PhiID_DC		= PhiID_Atoms.str + PhiID_Atoms.stx + ...
								PhiID_Atoms.sty;
							PhiID_CD		= PhiID_Atoms.sts;
							PhiID_DoubleRed	= PhiID_Atoms.rtr;
							PhiID_UC		= PhiID_Atoms.rts + PhiID_Atoms.xts + ...
								PhiID_Atoms.yts;
							PhiID_Syn		= PhiID_UC + PhiID_DC + PhiID_Atoms.sts;
							PhiID_Transfer	= PhiID_Atoms.xty + PhiID_Atoms.ytx;

							if strcmp(measure, 'phiid_measures_mmi');
								results(h,i,j).PhiID_CE_MMI(k,p)		= PhiID_CE;
								results(h,i,j).PhiID_DC_MMI(k,p)		= PhiID_DC;
								results(h,i,j).PhiID_CD_MMI(k,p)		= PhiID_CD;
								results(h,i,j).PhiID_DoubleRed_MMI(k,p)	= PhiID_DoubleRed;
								results(h,i,j).PhiID_UC_MMI(k,p)		= PhiID_UC;
								results(h,i,j).PhiID_Syn_MMI(k,p)		= PhiID_Syn;
								results(h,i,j).PhiID_Transfer_MMI(k,p)	= PhiID_Transfer;

								fprintf('\n\tPhiID_CE_MMI = %g\n', PhiID_CE);

							elseif strcmp(measure, 'phiid_measures_ccs');
								results(h,i,j).PhiID_CE_CCS(k,p)		= PhiID_CE;
								results(h,i,j).PhiID_DC_CCS(k,p)		= PhiID_DC;
								results(h,i,j).PhiID_CD_CCS(k,p)		= PhiID_CD;
								results(h,i,j).PhiID_DoubleRed_CCS(k,p)	= PhiID_DoubleRed;
								results(h,i,j).PhiID_UC_CCS(k,p)		= PhiID_UC;
								results(h,i,j).PhiID_Syn_CCS(k,p)		= PhiID_Syn;
								results(h,i,j).PhiID_Transfer_CCS(k,p)	= PhiID_Transfer;

								fprintf('\n\tPhiID_CE_CCS = %g\n', PhiID_CE);
							end

						%------------------------------------------------------------------------
						% INTEGRATED INFORMATION MEASURES
						%------------------------------------------------------------------------
						elseif strcmp(measure, 'integrated_info_measures')

							% compute time-lagged covariance matrix of transformed
							% KM time series
							phases_autocov_sequence = tsdata_to_autocov(...
								phases_transformed, 1);

							% construct full time-lagged covariance matrix
							% (block-toeplitz form)
							phases_lag01_autocov = [phases_autocov_sequence(:,:,1), ...
								phases_autocov_sequence(:,:,2)'; ...
								phases_autocov_sequence(:,:,2), ...
								phases_autocov_sequence(:,:,1)];

							for z = 1:length(list_integrated_info_measures);

								% name template to instantiate JIDT calculators
								class_template = ['infodynamics_with_iit.measures.continuous.' ...
									'gaussian.%sCalculatorGaussian'];

								% the following uses only the full time-lagged covariance
								% matrix to compute the measure
								if strcmp(list_integrated_info_measures{z}, ...
										'AverageCorrelation')
									diag_phases_lag01_autocov = diag(diag(...
										phases_lag01_autocov));

									corrmat = (diag_phases_lag01_autocov^-0.5) * ...
										diag_phases_lag01_autocov * ...
										(diag_phases_lag01_autocov^-0.5);

									measure_value = (sum(abs(corrmat(:))) - ...
										sum(diag(abs(corrmat))))/(n_nodes^2-n_nodes);
								else

									try
										calc = javaObject(sprintf(class_template, ...
											list_integrated_info_measures{z}));
										calc.initialise(n_nodes);
										calc.setLaggedCovariance(...
											phases_lag01_autocov);
										measure_value = ...
											calc.computeAverageLocalOfObservations();
									catch
										measure_value = NaN;
									end
								end

								switch list_integrated_info_measures{z}
									case 'AverageCorrelation',    results(i,j). ...
											AverageCorr(k,p)		   = ...
											measure_value;
									case 'IntegratedInformation', results(i,j). ...
											IntegratedInfo(k,p)	   = ...
											measure_value;
									case 'IntegratedInteraction', results(i,j). ...
											IntegratedInteraction(k,p) = ...
											measure_value;
									case 'DecoderIntegration',    results(i,j). ...
											DecoderIntegration(k,p)	   = ...
											measure_value;
									case 'CausalDensity',         results(i,j). ...
											CausalDensity(k,p)	   = ...
											measure_value;
									case 'IntegratedSynergy',     results(i,j). ...
											IntegratedSynergy(k,p)	   = ...
											measure_value;
									case 'TimeDelayedMutualInfo', results(i,j). ...
											TimeDelayedMI(k,p)         = ...
											measure_value;
									otherwise, error('Unknown ''calcName''');
								end % switch list_integrated_info_measures{z}

								fprintf('\n\tmeasure value = %g\n', measure_value);
							end


						end
					otherwise error('Unknown measure ''%s''', measure);
					end % switch(measure)

					clear noise_cov;
					clear phases_autocov_sequence;
					clear phases;
					clear phases_transformed;
					clear phases_rica
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

if strcmp(measure, 'phiid_measures_mmi') || strcmp(measure, ...
		'phiid_measures_ccs')

 	if strcmp(measure, 'phiid_measures_mmi')
		phiid_fieldnames = {'PhiID_DoubleRed_MMI', ...
			'PhiID_CE_MMI', 'PhiID_DC_MMI', ...
			'PhiID_CD_MMI', 'PhiID_UC_MMI', ...
			'PhiID_Syn_MMI', 'PhiID_Transfer_MMI'};

	elseif strcmp(measure, 'phiid_measures_ccs')
                phiid_fieldnames = {'PhiID_DoubleRed_CCS', ...
			    'PhiID_CE_CCS', 'PhiID_DC_CCS', ...
			    'PhiID_CD_CCS', 'PhiID_UC_CCS', ...
			    'PhiID_Syn_CCS', 'PhiID_Transfer_CCS'};
	end

	% must correspond to order of PhiID-based measures in
	% 'phiid_fieldnames'
	phiid_filenames = {'double_red', 'ce', 'dc' 'cd', ...
		'uc', 'syn', 'transfer'};

 	for h = 1:length(phiid_fieldnames)
		for i=1:size(results,1)
			for j=1:size(results,2)

			results_phiid(i,j).n_nodes		     = ...
				results(i,j).n_nodes;
			results_phiid(i,j).model		     = ...
				results(i,j).model;
			results_phiid(i,j).n_nodes		     = ...
				results(i,j).n_nodes;
			results_phiid(i,j).spectral_radius	     = ...
				results(i,j).spectral_radius;
			results_phiid(i,j).cross_coup_range	     = ...
				results(i,j).cross_coup_range;
			results_phiid(i,j).rmi_range		     = ...
				results(i,j).rmi_range;
			results_phiid(i,j).global_coupling_mag   = ...
				results(i,j).global_coupling_mag;
			results_phiid(i,j).noise_corr	           = ...
				results(i,j).noise_corr;
			results_phiid(i,j).coupling_matrices     = ...
				results(i,j).coupling_matrices;
			results_phiid(i,j).(phiid_fieldnames{h}) = ...
				results(i,j).(phiid_fieldnames{h});

			end
		end

		% save results
		results_phiid_filename = [results_filename '_' ...
			phiid_filenames{h}'];
		results_filepath = fullfile(pathout_data_measures, ...
			[model_specific_filename '_' ...
			results_phiid_filename '.mat']);
		fprintf('\nSaving results file ''%s'' ...', ...
			results_filepath);
		save(results_filepath, 'results_phiid', '-v7.3');
		fprintf(' one\n');

	end

else

	% save results
	results_filepath = fullfile(pathout_data_measures, ...
		[model_specific_filename '_' results_filename '.mat']);
	fprintf('\nSaving results file ''%s'' ...', ...
		results_filepath);
	save(results_filepath, 'results', '-v7.3');
	fprintf(' one\n');
end


