%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mec_random_cross_coup_rmi.m
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
% >> params_dd_ce_co_info_rand_cross_coup_rmi.m;
% >> mec_mvar_rand_cross_coup_rmi.m;
%
% To run in batch mode, See Bash script 'bash_mec_rand_cross_coups_rmi.sh'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DISPLAY PARAMETERS

fprintf('\nParameters:\n-----------\n\n');
fprintf('current directory: %s\n',		pwd);
fprintf('base directory: %s\n',		basedir);
fprintf('measure: %s\n',			measure);
fprintf('n_nodes: %d\n',			n_nodes);
fprintf('ar_model_order: %d\n',		ar_model_order);
fprintf('spectral_radius: %g\n',		spectral_radius);
fprintf('n_samples_cross_coup: %d\n',	n_samples_cross_coup);
fprintf('n_samples_noise_corr: %d\n',	n_samples_noise_corr);
fprintf('cross_coup_range: %g - %g\n',	cross_coup_range(1), ...
	cross_coup_range(2));
fprintf('n_cross_coup: %d\n',			n_cross_coup);
fprintf('rmi_range: %g - %g\n',		rmi_range(1), rmi_range(2));
fprintf('n_rmi: %d\n',				n_rmi);
fprintf('results_filename: %s\n',		results_filename);
fprintf('\n');

% only for DD/CE/CO-INFO
if exist('m_dim_ce_dd_co_info')
	fprintf('m_dim_ce_dd_co_info: %d\n', m_dim_ce_dd_co_info);
end
if exist('dim_reduction')
	fprintf('dim_reduction: %s\n', strjoin(dim_reduction,', '));
end
if exist('n_samples_dd_ce_co_info')
	fprintf('n_samples_dd_ce_co_info: %d\n', n_samples_dd_ce_co_info);
end

% only for PhiID_CE_CCS and other PhiID-based measures using CCS
if exist('n_datapoints_phiid_ccs')
	fprintf('n_datapoints_phiid_ccs: %d\n', n_datapoints_phiid_ccs);
end

%% SET PATHS

% run SSDI startup (also calls MVGC2 startup)
run(fullfile(basedir, 'packages_and_code_repos/', 'ssdi', 'startup'));

addpath(fullfile(basedir, 'functions'));
addpath(fullfile(basedir, 'scripts', 'multivar_autoreg', ...
	'rand_cross_coup_rmi'));

% {
javaaddpath(fullfile(basedir, 'scripts', 'multivar_autoreg', ...
	'rand_cross_coup_rmi', 'infodynamics_with_iit.jar'));
javaaddpath(fullfile(basedir, 'scripts', 'multivar_autoreg', ...
	'rand_cross_coup_rmi', 'commons-math3-3.5.jar'));
javaaddpath(fullfile(basedir, 'functions', 'infodynamics_with_iit.jar'));
javaaddpath(fullfile(basedir, 'functions', 'commons-math3-3.5.jar'));

try
	javaaddpath(fullfile(basedir, 'scripts', 'multivar_autoreg', ...
		'rand_cross_coup_rmi', 'infodynamics_with_iit.jar'));
	disp('Path added successfully.');
catch exception
	disp(['Error adding path: ' exception.message]);
end

javaClassPath = javaclasspath;
disp(javaClassPath);
%}

rng(seed) % seed random number generator

% loop over coupling matrices and RMI values (including random sampling of 
% connectivities & noise correlations)

%% MEMORY PRE-ALLOCATION 

% use struct array for results

% coupling/noise zero matrices
coup_noise_corr_zeros = zeros(n_samples_cross_coup, n_samples_noise_corr);
noise_corr_zeros = zeros(n_nodes, n_nodes, n_samples_cross_coup, ...
        n_samples_noise_corr);
coup_zeros = zeros(n_nodes, n_nodes, ar_model_order, ...
        n_samples_cross_coup, n_samples_noise_corr);

% the same for all measures
blank_struct.model			= model_specific_filename;
blank_struct.n_nodes			= n_nodes;
blank_struct.spectral_radius		= spectral_radius;
blank_struct.global_coup		= coup_noise_corr_zeros;
blank_struct.noise_corr			= noise_corr_zeros;
blank_struct.coupling_matrices		= coup_zeros;
blank_struct.rmi_for_this_sim		= 0;
blank_struct.cross_coup_for_this_sim	= 0;

% cross-coupling (exponential scale)
all_cross_coup = 2.^linspace(cross_coup_range(1), ...
	cross_coup_range(2), n_cross_coup); 
% rmi values  (exponential scale)
all_rmi = 2.^linspace(rmi_range(1), rmi_range(2), n_rmi); 

blank_struct.cross_coup_range		= [all_cross_coup(1), all_cross_coup(end)];
blank_struct.rmi_range			= [all_rmi(1), all_rmi(end)];

% measure-dependent
switch measure

	case 'multi_info'
		blank_struct.MultiInfo = coup_noise_corr_zeros;

	case 'dd_ce_co_info'
		blank_struct.m_dim_ce_dd_co_info	= m_dim_ce_dd_co_info;

		if any(strcmp(dim_reduction,'pca'))
			blank_struct.DD_PCA		= coup_noise_corr_zeros;
			blank_struct.ShannonCE_PCA	= coup_noise_corr_zeros;
			blank_struct.CoInfoPCA		= coup_noise_corr_zeros;
		end

		if any(strcmp(dim_reduction,'grassmanian'))
			blank_struct.DDGrassMin			= coup_noise_corr_zeros;
			blank_struct.DDGrassMean		= coup_noise_corr_zeros;
			blank_struct.ShannonCEGrassMax		= coup_noise_corr_zeros;
			blank_struct.ShannonCEGrassMean		= coup_noise_corr_zeros;
			blank_struct.CoInfoGrassMin		= coup_noise_corr_zeros;
			blank_struct.CoInfoGrassMax		= coup_noise_corr_zeros;
			blank_struct.CoInfoGrassMean		= coup_noise_corr_zeros;
		end

	case 'phiid_measures_mmi'
		blank_struct.PhiID_DoubleRed_MMI	= coup_noise_corr_zeros;
		blank_struct.PhiID_CE_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_DC_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_CD_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_UC_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_Syn_MMI		= coup_noise_corr_zeros;
		blank_struct.PhiID_Transfer_MMI		= coup_noise_corr_zeros;

	case 'phiid_measures_ccs'
		blank_struct.PhiID_DoubleRed_CCS	= coup_noise_corr_zeros;
		blank_struct.PhiID_CE_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_DC_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_CD_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_UC_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_Syn_CCS		= coup_noise_corr_zeros;
		blank_struct.PhiID_Transfer_CCS		= coup_noise_corr_zeros;

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
	
	case 'control_measures'
		if any(strcmp(list_control_measures, 'AverageCorrelation'   )), ...
				blank_struct.AverageCorr		= coup_noise_corr_zeros; end
		if any(strcmp(list_control_measures, 'TimeDelayedMutualInfo')), ...
				blank_struct.TimeDelayedMI		= coup_noise_corr_zeros; end

	otherwise error('Unknown measure ''%s''', measure);
end % switch measure
	
results = repmat(blank_struct, n_cross_coup, n_rmi); % preallocate

w = whos('results');
fprintf('Pre-allocated size of ''results'' = %.4f MB\n', ...
	w.bytes/1000/1000);

%% SET UP MVAR MODEL 

% logical indices of on/off-diagonal matrix elements for VAR coefficients 
% sequence
var_coeff_on_diag = repmat(diag(true(n_nodes, 1)), 1, 1, ...
	ar_model_order);
var_coeff_off_diag = ~var_coeff_on_diag;
i
%% MODEL PARAMETER SWEEP (RMI & CROSS-COUPLING)

for i = 1:n_cross_coup;
	fprintf('\nCROSS-COUPLING %d of %d = %g\n', i, n_cross_coup, ...
		all_cross_coup(i));

	for j = 1:n_rmi;
		fprintf('\n\tRMI %d of %d = %g\n', j, n_rmi, all_rmi(j));
		fprintf('\nCROSS-COUPLING %d of %d = %g\n', i, n_cross_coup, ...
			all_cross_coup(i));
		st = tic;

		results(i,j).cross_coup_this_sim	= all_cross_coup(i);
		results(i,j).rmi_for_this_sim		= all_rmi(j);

		for k = 1:n_samples_cross_coup;
			fprintf('\t\tcross-couplings sample %d of %d\n', ...
				k, n_samples_cross_coup);
				
			for p = 1:n_samples_noise_corr;
				fprintf('\t\t\tnoise corr sample %d of %d\n', ...
					p, n_samples_noise_corr);

				%------------------------------------------------------------------------
				% COUPLING & NOISE CORRELATION MATRICES
				%------------------------------------------------------------------------
	
				var_coeff_seq = var_rand(n_nodes, ar_model_order, 1);

				% scale all off-diagonal elements
				var_coeff_seq(var_coeff_off_diag) = all_cross_coup(i)* ...
					var_coeff_seq(var_coeff_off_diag);

				% enforce spectral norm as specified
				var_coeff_seq = specnorm(var_coeff_seq, spectral_radius);

				results(i,j).global_coup(k,p) = sqrt(mean(...
					var_coeff_seq(var_coeff_off_diag).^2)/ ...
					mean(var_coeff_seq(:).^2));

				% generate random noise correlation matrices
				% note: -rmi sets rmi relative to uniform (i.e., size-adjusted) 
				% correlation matrix at the given dimension
				noise_corr = corr_rand(n_nodes, -all_rmi(j));

				results(i,j).noise_corr(:,:,k,p) = noise_corr;
				results(i,j).coupling_matrices(:,:,:,k,p) = var_coeff_seq;

				% calculate autocovariance sequence (up to 1 lag only)
				% -1 means *exactly* one lag!
				autocov_seq = var_to_autocov(var_coeff_seq, noise_corr, -1); 

				covmat  = autocov_seq(:,:,1);			% unlagged covariance matrix
				lag1_autocov = autocov_seq(:,:,2);		% 1-lag autocovariance
				lag01_autocov = [covmat, lag1_autocov; ...	% 0,1-lag autocovariance 
					lag1_autocov', covmat];			% (block-Toeplitz matrix)

				switch measure

					%------------------------------------------------------------------------
					% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
					%------------------------------------------------------------------------
					case 'dd_ce_co_info'

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
								results(i,j).DDGrassMin(k,p)		= min(DD);
								results(i,j).DDGrassMean(k,p)		= mean(DD);
								results(i,j).ShannonCEGrassMax(k,p)	= max(CE);
								results(i,j).ShannonCEGrassMean(k,p)	= mean(CE);
								results(i,j).CoInfoGrassMin(k,p)	= min(CI);
								results(i,j).CoInfoGrassMax(k,p)	= max(CI);
								results(i,j).CoInfoGrassMean(k,p)	= mean(CI);

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

						% need to compute local entropies for CCS, hence simulating
                                    		% MVAR time-series
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
							class_template = ['infodynamics_with_iit.measures.continuous.' ...
								'gaussian.%sCalculatorGaussian'];

							% the following uses only the full time-lagged covariance 
							% matrix to compute the measure
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
						

							switch list_integrated_info_measures{z}

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
								otherwise, error('Unknown integrated information measure.');

							end % switch list_integrated_info_measures{z}
							fprintf('\n\tmeasure value = %g\n', measure_value);
						end
					
					%------------------------------------------------------------------------
					% CONTROL MEASURES
					%------------------------------------------------------------------------
					case 'control_measures'

						for z = 1:length(list_control_measures);

							if strcmp(list_control_measures{z}, ...
									'AverageCorrelation')
								%diag_lag01_autocov = diag(diag(lag01_autocov));		
								%corrmat = (diag_lag01_autocov^-0.5) * ...
								%	diag_lag01_autocov * (diag_lag01_autocov^-0.5);
								%measure_value = (sum(abs(corrmat(:))) - ...
                                                                %        sum(diag(abs(corrmat))))/(n_nodes^2-n_nodes);

								% average covariance matrix
								corrmat = corrcov(covmat);
								average_corr = mean(tril(corrmat,-1), 'all');

								results(i,j).AverageCorr(k,p) = average_corr;
								fprintf('\n\tAverage Corr = %g\n', average_corr);

							elseif strcmp(list_control_measures{z}, ...
									'TimeDelayedMutualInfo')
								
								% name template to instantiate JIDT calculators
                                                        	class_template = ['infodynamics_with_iit.measures.continuous.' ...
                                                                'gaussian.%sCalculatorGaussian'];

								% the following uses only the full time-lagged covariance
                                                        	% matrix to compute the measure

								%calc = javaObject(sprintf(class_template, ...
                                                                %        list_control_measures{z}))
                                                                %calc.initialise(n_nodes)
                                                                %calc.setLaggedCovariance(lag01_autocov)
                                                                %tdmi = calc.computeAverageLocalOfObservations()

								try
									calc = javaObject(sprintf(class_template, ...
										list_control_measures{z}));
									calc.initialise(n_nodes); 
									calc.setLaggedCovariance(lag01_autocov);

									tdmi = calc.computeAverageLocalOfObservations();
									
								catch
									tdmi = NaN;
								end
								results(i,j).TimeDelayedMI(k,p) = tdmi;
								fprintf('\n\tTDMI = %g\n', tdmi);
							
							elseif strcmp(list_control_measures{z}, 'MultiInfo')

                                                		% process multi-information (nats)
                                                		results(i,j).MultiInfo(k,p) = sum(log(diag(covmat))) - ...
                                                        		logdet(covmat);

                                                		fprintf('\n\tMulti Info = %g\n', results(i,j).MultiInfo(k,p));
							else 
								error('Unknown control measure.');
							end
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
			
			results_phiid(i,j).n_nodes		     	= ...
				results(i,j).n_nodes;
			results_phiid(i,j).model		     	= ...
				results(i,j).model;
			results_phiid(i,j).n_nodes		     	= ...
				results(i,j).n_nodes;
			results_phiid(i,j).spectral_radius	     	= ...
				results(i,j).spectral_radius;
			results_phiid(i,j).cross_coup_range	     	= ...
				results(i,j).cross_coup_range;
			results_phiid(i,j).rmi_range			= ...
				results(i,j).rmi_range;
			results_phiid(i,j).global_coup   		= ...
				results(i,j).global_coup;
			results_phiid(i,j).noise_corr	        	= ...
				results(i,j).noise_corr;
			results_phiid(i,j).coupling_matrices     	= ...
				results(i,j).coupling_matrices;
			results_phiid(i,j).(phiid_fieldnames{h}) 	= ...
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

