%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mec_metacommkm_rand_phase_lag_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this script, set 'basedir' and 'measure', and run the appropriate 
% parameter script
%
% E.g., to run interactively:
%
% >> basedir = '/its/home/ns508';
% >> measure = 'dd_ce_co_info';
% >> params_dd_ce_co_info_rand_phase_lag_rmi.m;
% >> mec_metacommkm_rand_phase_lag_rmi.m;
%
% To run in batch mode, See Bash script 'bash_mec_rand_cross_coups_rmi.sh'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DISPLAY PARAMETERS

fprintf('\nParameters:\n-----------\n\n');
fprintf('current directory: %s\n',			pwd);
fprintf('base directory: %s\n',			basedir);
fprintf('results_filename: %s\n',			results_filename);
fprintf('measure: %s\n',				measure);
fprintf('ar_model_order: %d\n',			ar_model_order);
fprintf('seed: %d\n',					seed);
% -------------------------------------------------------------------------
fprintf('n_nodes: %d\n',				n_nodes);
fprintf('n_comms: %d\n',				n_comms);
fprintf('n_inter_comm_coups: %d\n',			n_inter_comm_coups);
fprintf('n_nodes_per_comm: %d\n',			n_nodes_per_comm);
% -------------------------------------------------------------------------
fprintf('phase_lag_range: %g - %g\n',		phase_lag_range(1), ...
	phase_lag_range(2));
fprintf('n_phase_lag: %d\n',				n_phase_lag);
% -------------------------------------------------------------------------
fprintf('rmi_range: %g - %g\n',			rmi_range(1), rmi_range(2));
fprintf('n_rmi: %d\n',					n_rmi);
fprintf('n_samples_noise_corr: %d\n',		n_samples_noise_corr);
fprintf('noise_mag: %d\n',				noise_mag);
fprintf('adjust_rmi: %d\n',				adjust_rmi);
% -------------------------------------------------------------------------
fprintf('coup_ratio_mean_range: %g - %g\n',	coup_ratio_mean_range(1), ...
	coup_ratio_mean_range(2));
fprintf('n_coup_ratio_mean: %d\n',			n_coup_ratio_mean);
fprintf('n_samples_coup: %d\n',			n_samples_coup);
fprintf('coup_ratio_mean_rel_dev: %d\n',		coup_ratio_mean_rel_dev);
fprintf('global_coup: %d\n',				global_coup);
% -------------------------------------------------------------------------
fprintf('n_sim_points: %d\n',				n_sim_points);
fprintf('n_equil_points: %d\n',			n_equil_points);
fprintf('dt_sim_points: %d\n',			dt_sim_points);
fprintf('sim_mode: %d\n',				sim_mode);
% -------------------------------------------------------------------------
fprintf('n_local_gc: %d\n',				n_local_gc);
fprintf('emp_sample: %d\n',				emp_sample);
fprintf('stabilise_care: %d\n',			stabilise_care);
% -------------------------------------------------------------------------
fprintf('osc_freq_mean: %d\n',			osc_freq_mean);
fprintf('osc_freq_rel_dev: %d\n',			osc_freq_rel_dev);
% -------------------------------------------------------------------------
fprintf('\n');

% --------------------------------------------
% measure-specific parameters
% --------------------------------------------

% only for DD/CE/CO-INFO
if exist('n_samples_dd_ce_co_info')
	fprintf('n_samples_dd_ce_co_info: %d\n', n_samples_dd_ce_co_info);
end

%% SET PATHS

% run kvar (also calls SSDI & MVGC2 startup)
run(fullfile(basedir,'packages_and_code_repos', 'kvar', 'startup'));

addpath(fullfile(basedir, 'code', 'common'));
addpath(fullfile(basedir, 'code', 'kuramoto_osc', ...
	'meta_comm_rand_phase_lag_rmi', 'functions'));
addpath(fullfile(basedir, 'code', 'kuramoto_osc', ...
	'meta_comm_rand_phase_lag_rmi', 'scripts'));

% {
javaaddpath(fullfile(basedir, '/infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'kuramoto_osc', ...
	'meta_comm_rand_phase_lag_rmi', 'scripts', 'infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'kuramoto_osc', ...
	'meta_comm_rand_phase_lag_rmi', 'functions', 'infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'kuramoto_osc', ...
	'meta_comm_rand_phase_lag_rmi', 'scripts', 'commons-math3-3.5.jar'));
javaaddpath(fullfile(basedir, 'code', 'common', 'infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'common', 'commons-math3-3.5.jar'));

try
	javaaddpath(fullfile(basedir, 'code', 'kuramoto_osc', ...
	'meta_comm_rand_phase_lag_rmi', 'scripts', 'infodynamics.jar'));
	disp('Path added successfully.');
catch exception
	disp(['Error adding path: ' exception.message]);
end

javaClassPath = javaclasspath;
disp(javaClassPath);
%}

rng(seed) % seed random number generator

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
blank_struct.n_inter_comm_coups		= n_inter_comm_coups;
blank_struct.n_nodes_per_comm			= n_nodes_per_comm;
blank_struct.ar_model_order			= ar_model_order;

blank_struct.phase_lag_range			= phase_lag_range;

blank_struct.rmi_range				= rmi_range;
blank_struct.noise_mag				= noise_mag;
blank_struct.noise_corr				= noise_corr_zeros;

blank_struct.coup_ratio_mean_range		= coup_ratio_mean_range;
blank_struct.coup_ratio_mean_rel_dev	= coup_ratio_mean_rel_dev;
blank_struct.global_coup			= global_coup;
blank_struct.coupling_matrices		= coup_zeros;

blank_struct.n_sim_points			= n_sim_points;
blank_struct.n_equil_points			= n_equil_points;
blank_struct.dt_sim_points			= dt_sim_points;
blank_struct.sim_mode				= sim_mode;

blank_struct.n_local_gc				= n_local_gc;
blank_struct.emp_sample				= emp_sample;
blank_struct.stabilise_care			= stabilise_care;

blank_struct.osc_freq_mean			= osc_freq_mean;
blank_struct.osc_freq_rel_dev			= osc_freq_rel_dev;

blank_struct.seed	= seed;

% measure-dependent
switch measure

	case 'multi_info'
		blank_struct.MultiInfo = coup_noise_corr_zeros;

	case 'dd_ce_co_info'
		if any(strcmp(dim_reduction,'comm_sync'))
			blank_struct.DDCommSync			= coup_noise_corr_zeros;
			blank_struct.ShannonCECommSync	= coup_noise_corr_zeros;
			blank_struct.CoInfoCommSync		= coup_noise_corr_zeros;
			blank_struct.m_dim_ce_dd_co_info	= n_comms;
		else 
			error('Undefined dim_reduction, must be "comm_sync".')
		end

	case 'phiid_measures_mmi'
        blank_struct.PhiID_DoubleRed_MMI	= coup_noise_corr_zeros;
        blank_struct.PhiID_CE_MMI		= coup_noise_corr_zeros;
        blank_struct.PhiID_DC_MMI		= coup_noise_corr_zeros;
        blank_struct.PhiID_CD_MMI		= coup_noise_corr_zeros;
        blank_struct.PhiID_UC_MMI		= coup_noise_corr_zeros;
        blank_struct.PhiID_Syn_MMI		= coup_noise_corr_zeros;
        blank_struct.PhiID_Transfer_MMI	= coup_noise_corr_zeros;

	case 'phiid_measures_ccs'
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
				blank_struct.DecoderIntegration	= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'CausalDensity'        )), ...
				blank_struct.CausalDensity		= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'IntegratedSynergy'    )), ...
				blank_struct.IntegratedSynergy	= coup_noise_corr_zeros; end
		if any(strcmp(list_integrated_info_measures, ...
				'TimeDelayedMutualInfo')), ...
				blank_struct.TimeDelayedMI		= coup_noise_corr_zeros; end

	otherwise error('Unknown measure ''%s''', measure);
end % switch measure

%% CONTINUE HERE
%--------------------------------------------------------------------------
% CONTINUE HERE
% How do we need to change 'results' to accommodate3 variables?
%--------------------------------------------------------------------------

results = repmat(blank_struct, n_coup_ratio_mean, n_rmi); % preallocate

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
all_phase_lag = pi/2 - linspace(phase_lag_range(1), ...		% phase lag
	phase_lag_range(2), n_phase_lag); 
all_coup_ratio_mean = linspace(coup_ratio_mean_range(1), ...	% coupling ratio mean values
	coup_ratio_mean_range(2), n_coup_ratio_mean);

n_nodes = n_comms * n_nodes_per_comm;					% total number of oscillators
osc_freq_dev = osc_freq_rel_dev * osc_freq_mean;			% oscillator frequencies deviation
coup_ratio_mean_dev = coup_ratio_mean_rel_dev * coup_ratio_mean;	% coupling ratio mean deviation

% oscillator natural frequencies
if osc_freq_dev > 0
	seed_state = rng_seed(seed_osc_freq);
	osc_freq = pi*ms_betarnd(osc_freq_mean/pi, ...			% oscillator frequencies Beta 
		osc_freq_dev, n_nodes, 1);					% distributed on [0,pi] with  
	rng_restore(seed_state);						% mean osc_freq_mean and 
else											% std. dev osc_freq_dev
	osc_freq = osc_freq_mean * ones(n_nodes, 1);
end
sampling_freq  = 1/dt_sim_points;						% sampling frequency
nyquist_freq = sampling_freq/2;						% Nyqvist frequency
freq_hz = nyquist_freq * osc_freq/pi;					% frequencies in Hz

% simulation times
dt_equil_points = floor(n_equil_points / dt_sim_points);		% equilibriation time increments
assert(dt_equil_points >= 0, 'Equilibriation time too short,' ...
	' or time increment too large!');
n_equil_points = dt_equil_points * dt_sim_points;			% adjusted equilibriation time

sim_time_increments = floor(n_sim_points/dt_sim_points);		% simulation time increments
assert(sim_time_increments > 0,'simulation time too short,' ...
	' or time increment too large!');
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
	if coup_ratio_mean_rel_dev = 0 && n_samples_coup = 1
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
			fprintf('\nCOUPLING %d of %d = %g\n', h, n_coup_ratio_mean, ...
				all_coup_ratio_mean(h));
			fprintf('\nPHASE LAG %d of %d = %g\n', i, n_phase_lag, ...
				all_phase_lag(i));
			fprintf('\n\tRMI %d of %d = %g\n', j, n_rmi, all_rmi(j));
			st = tic;

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
				assert(~isempty(coupling_matrix), 'Coupling matrix does not exist.');

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
					noise_corr = randn(tot_sim_time_increments, n_nodes) * chol(noise_cov);
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
						noise_corr, ...
						sim_mode);

					% truncate equilibriation
					phases = phases(:, dt_equil_points+1:end);
					order_param = order_param(dt_equil_points+1:end);

					% community order parameters (synchronisation)
					% --> macroscopic variable
					order_param_comms = zeros(n_comms, sim_time_increments);
					for k = 1:n_comms
						order_param_comms(k,:) = hypot(mean(cos(phases(communities{k},:))), ...
							mean(sin(phases(communities{k},:))));
					end

					% ---------------------------------
					% AMEND TO ACCOMMODATE THREE PARAMS
					% ---------------------------------
					results(h,i,j).noise_corr(:,:,k,p) = noise_corr;
					results(h,i,j).coupling_matrices(:,:,:,k,p) = coupling_matrix;

					results(h,i,j).global_coup(:,:,:,k,p) = sqrt(mean(...
						coupling_matrix(diag(coupling_matrix)).^2)/ ...
						mean(coupling_matrix(:).^2));

				switch measure

					%------------------------------------------------------------------------
					% MULTI-INFORMATION
					%------------------------------------------------------------------------
					case 'multi_info'
		
						% process multi-information (nats)
						results(i,j).MultiInfo(k,p) = sum(log(diag(covmat))) - ...
							logdet(covmat); 

						fprintf('\n\tMultiInfo = %g\n', results(i,j).MultiInfo(k,p));

					%------------------------------------------------------------------------
					% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
					%------------------------------------------------------------------------
					case 'dd_ce_co_info'
						
						% -------------------
						% TODO
						m_dim_ce_dd_co_info     = 1;
						dim_reduction         	= {'pca','grassmanian'};
						% -------------------
						
						% convert from VAR to ISS
						[state_trans_mat, observation_mat, kalman_gain_mat] = var_to_ss( ...
							var_coeff_seq, noise_corr);

						% pre-compute state prediction error covariance matrices
						[~, error_cov] = iss2ce_precomp(state_trans_mat, observation_mat, ...
							kalman_gain_mat, noise_corr);

						% calculate DD & CE
						for z = 1:length(dim_reduction)

							if strcmp(dim_reduction{z}, 'pca')

								[~, macro]	= pcacov(covmat);
								[CI, DD]	= iss2ce(macro, state_trans_mat, ...
									observation_mat, kalman_gain_mat, noise_corr, ...
									covmat, error_cov);
								CE		= CI - DD;

								results(i,j).DD_PCA(k,p)		= DD;
								results(i,j).ShannonCE_PCA(k,p)	= CE;
								results(i,j).CoInfoPCA(k,p)		= CI;

								fprintf('\n\tCoInfoPCA = %g\n', results(i,j). ...
									CoInfoPCA(k,p));

							elseif strcmp(dim_reduction{z}, 'grassmanian')

								macro	= rand_orthonormal(n_nodes, ...
									m_dim_ce_dd_co_info, n_samples_dd_ce_co_info);
								CI	= zeros(n_samples_dd_ce_co_info,1);
								DD	= zeros(n_samples_dd_ce_co_info,1);

								for s = 1:n_samples_dd_ce_co_info
									[CI(s),DD(s)] = iss2ce(macro(:,:,s), ...
										state_trans_mat, observation_mat, ...
										kalman_gain_mat, noise_corr, covmat, ...
										error_cov);
								end

								CE = CI - DD;
								results(i,j).DDGrassMin(k,p)			= min(DD);
								results(i,j).DDGrassMean(k,p)			= mean(DD);
								results(i,j).ShannonCEGrassMax(k,p)		= max(CE);
								results(i,j).ShannonCEGrassMean(k,p)	= mean(CE);
								results(i,j).CoInfoGrassMin(k,p)		= min(CI);
								results(i,j).CoInfoGrassMax(k,p)		= max(CI);
								results(i,j).CoInfoGrassMean(k,p)		= mean(CI);

								fprintf('\n\tCoInfoGrassMax = %g\n', ...
									results(i,j).CoInfoGrassMax(k,p));

							end

						end
					
					%------------------------------------------------------------------------
                              % PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, 
					% CAUSAL DECOUPLING USING MMI
                              %------------------------------------------------------------------------
					case 'phiid_measures_mmi'

						PhiID_Atoms		= PhiIDFullAnalytical(lag01_autocov, 'MMI');

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

						%results_phiid_atoms(i,j).PhiID_Atoms_MMI(k,p) = PhiID_Atoms_MMI;
						results(i,j).PhiID_CE_MMI(k,p)		= PhiID_CE;
						results(i,j).PhiID_DC_MMI(k,p)		= PhiID_DC;
						results(i,j).PhiID_CD_MMI(k,p)		= PhiID_CD;
						results(i,j).PhiID_DoubleRed_MMI(k,p)	= PhiID_DoubleRed;
						results(i,j).PhiID_UC_MMI(k,p)		= PhiID_UC;
						results(i,j).PhiID_Syn_MMI(k,p)		= PhiID_Syn;
						results(i,j).PhiID_Transfer_MMI(k,p)	= PhiID_Transfer;

						fprintf('\n\tPhiID_CE_MMI = %g\n', PhiID_CE);

					%------------------------------------------------------------------------
                              % PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, 
					% CAUSAL DECOUPLING USING CCS
                              %------------------------------------------------------------------------						    
					case 'phiid_measures_ccs'

						[micro, corr_micro] = sim_mvar_model_with_lagged_cov(...
							n_datapoints_phiid_ccs, lag01_autocov, var_coeff_seq, ...
							ar_model_order);

						PhiID_Atoms		= PhiIDFullContinuousData(micro, ...
							ar_model_order, 'CCS', 'Gaussian');
						%PhiID_Atoms		= PhiIDFull_Analytical(lag01_autocov, 'CCS');

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

						%results_phiid_atoms(i,j).PhiID_Atoms_MMI(k,p) = PhiID_Atoms_MMI;
						results(i,j).PhiID_CE_CCS(k,p)		= PhiID_CE;
						results(i,j).PhiID_DC_CCS(k,p)		= PhiID_DC;
						results(i,j).PhiID_CD_CCS(k,p)		= PhiID_CD;
						results(i,j).PhiID_DoubleRed_CCS(k,p)	= PhiID_DoubleRed;
						results(i,j).PhiID_UC_CCS(k,p)		= PhiID_UC;
						results(i,j).PhiID_Syn_CCS(k,p)		= PhiID_Syn;
						results(i,j).PhiID_Transfer_CCS(k,p)	= PhiID_Transfer;

						fprintf('\n\tPhiID_CE_CCS = %g\n', PhiID_CE);

					%------------------------------------------------------------------------
					% INTEGRATED INFORMATION MEASURES
					%------------------------------------------------------------------------
					case 'integrated_info_measures'

						for z = 1:length(list_integrated_info_measures);

							% name template to instantiate JIDT calculators
							class_template = ['infodynamics.measures.continuous.' ...
								'gaussian.%sCalculatorGaussian'];

							% the following uses only the full time-lagged covariance 
							% matrix to compute the measure
							if strcmp(list_integrated_info_measures{z}, ...
									'AverageCorrelation')
								diag_lag01_autocov = diag(diag(lag01_autocov));

								corrmat = (diag_lag01_autocov^-0.5) * ...
									diag_lag01_autocov * (diag_lag01_autocov^-0.5);

								measure_value = (sum(abs(corrmat(:))) - ...
									sum(diag(abs(corrmat))))/(n_nodes^2-n_nodes);
							else

								try
									calc = javaObject(sprintf(class_template, ...
										list_integrated_info_measures{z}));
									calc.initialise(n_nodes); 
									calc.setLaggedCovariance(lag01_autocov);
									measure_value = ...
										calc.computeAverageLocalOfObservations();
								catch
									measure_value = NaN;
								end
							end

							switch list_integrated_info_measures{z}
								case 'AverageCorrelation',    results(i,j). ...
										AverageCorr(k,p)		   = measure_value;
								case 'IntegratedInformation', results(i,j). ...
										IntegratedInfo(k,p)	   = measure_value;
								case 'IntegratedInteraction', results(i,j). ...
										IntegratedInteraction(k,p) = measure_value;
								case 'DecoderIntegration',    results(i,j). ...
										DecoderIntegration(k,p)	   = measure_value;
								case 'CausalDensity',         results(i,j). ...
										CausalDensity(k,p)	   = measure_value;
								case 'IntegratedSynergy',     results(i,j). ...
										IntegratedSynergy(k,p)	   = measure_value;
								case 'TimeDelayedMutualInfo', results(i,j). ...
										TimeDelayedMI(k,p)         = measure_value;
								otherwise, error('Unknown ''calcName''');
							end % switch list_integrated_info_measures{z}

							fprintf('\n\tmeasure value = %g\n', measure_value);
						end

				otherwise error('Unknown measure ''%s''', measure);
				end % switch(measure)
			end
		end

		w = whos('results');
		et = toc(st);
		fprintf('\t\tcoupling/noise loop completed in %g seconds; size of ''results'' = %.4f MB\n', ...
			et, w.bytes/1000/1000);

	end
end

% save results
pathout_data_measures  = fullfile(basedir, 'results', 'analyses', ...
	'multivar_autoreg', 'rand_cross_coup_rmi', model_specific_path);

% row & column values for table
all_coup_ratio_mean_str = {};
for t = 1:length(all_coup_ratio_mean)
        all_coup_ratio_mean_str{t} = num2str(all_coup_ratio_mean(t));
end

all_rmi_str = {};
for e = 1:length(all_rmi)
        all_rmi_str{e} = num2str(all_rmi(e));
end


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
		
		% as table
		results_table_filename = [results_filename '_' ...
			phiid_filenames{h}'];
		results_table = array2table(results_phiid, 'RowNames', ...
			all_coup_ratio_mean_str, 'VariableNames', all_rmi_str);
		results_table_filepath = fullfile(pathout_data_measures, ...
			[model_specific_filename '_' ...
			results_table_filename '.mat']);
		fprintf('\nSaving results table file ''%s'' ...', ...
			results_table_filepath);
		save(results_table_filepath, 'results_table', '-v7.3');
		fprintf(' d1.0e+07 one\n');
		
	end

else 

	% save as table
	results_table_filename = results_filename;
	results_table = array2table(results, 'RowNames', ...
		all_coup_ratio_mean_str, 'VariableNames', all_rmi_str);
	results_table_filepath = fullfile(pathout_data_measures, ...
		[model_specific_filename '_' results_table_filename '.mat']);
	fprintf('\nSaving results table file ''%s'' ...', ...
		results_table_filepath);
	save(results_table_filepath, 'results_table', '-v7.3');
	fprintf(' d1.0e+07 one\n');
end


