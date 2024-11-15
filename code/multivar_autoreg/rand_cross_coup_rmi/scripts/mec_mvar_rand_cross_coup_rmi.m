%%%%%%%%%%%%%%%%1.0e+07 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mec_random_cross_coup_rmi.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this script, set 'basedir' and 'measure', and run the appropriate parameter script
%
% E.g., to run interactively:
%
% >> basedir = '/its/home/ns508';
% >> measure = 'dd_ce_co_info';
% >> mrc_parameters;
% >> mec_random_cross_coup_rmi;
%
% To run in batch mode, See Bash script 'bash_mec_rand_cross_coups_rmi.sh'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display parameters
fprintf('\nParameters:\n-----------\n\n');
fprintf('current directory: %s\n', pwd);
fprintf('base directory: %s\n', basedir);
fprintf('measure: %s\n', measure);
fprintf('n_nodes: %d\n', n_nodes);
fprintf('ar_model_order: %d\n', ar_model_order);
fprintf('spectral_radius: %g\n', spectral_radius);
fprintf('n_samples_cross_coup: %d\n', n_samples_cross_coup);
fprintf('n_samples_noise_corr: %d\n', n_samples_noise_corr);
fprintf('cross_coup_range: %g - %g\n',cross_coup_range(1), cross_coup_range(2));
fprintf('n_cross_coup: %d\n', n_cross_coup);
fprintf('rmi_range: %g - %g\n', rmi_range(1), rmi_range(2));
fprintf('n_rmi: %d\n', n_rmi);
fprintf('results_filename: %s\n', results_filename);
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

% run SSDI startup (also calls MVGC2 startup)
run(fullfile(basedir,'packages_and_code_repos', 'ssdi', 'startup'));

addpath(fullfile(basedir, 'packages_and_code_repos', 'ReconcilingEmergences-master'));
addpath(fullfile(basedir, 'code', 'common'));
addpath(fullfile(basedir, 'code', 'multivar_autoreg', 'rand_cross_coup_rmi', 'functions'));
addpath(fullfile(basedir, 'code', 'multivar_autoreg', 'rand_cross_coup_rmi', 'scripts'));

% {
javaaddpath(fullfile(basedir, '/infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'multivar_autoreg', 'rand_cross_coup_rmi', 'scripts', 'infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'multivar_autoreg', 'rand_cross_coup_rmi', 'functions', 'infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'multivar_autoreg', 'rand_cross_coup_rmi', 'scripts', 'commons-math3-3.5.jar'));
javaaddpath(fullfile(basedir, 'code', 'common', 'infodynamics.jar'));
javaaddpath(fullfile(basedir, 'code', 'common', 'commons-math3-3.5.jar'));

try
	javaaddpath(fullfile(basedir, 'code', 'multivar_autoreg', 'rand_cross_coup_rmi', 'scripts', 'infodynamics.jar'));
	disp('Path added successfully.');
catch exception
	disp(['Error adding path: ' exception.message]);
end

javaClassPath = javaclasspath;
disp(javaClassPath);
%}

rng(seed) % seed random number generator

% loop over coupling matrices and RMI values (including random sampling of connectivities & noise correlations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% memory pre-allocation (use struct array for results)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coupling/noise zero matrices
cnzeros	= zeros(n_samples_cross_coup, n_samples_noise_corr);
nncnzeros	= zeros(n_nodes, n_nodes, n_samples_cross_coup, ...
	n_samples_noise_corr);
nnmcnzeros	= zeros(n_nodes, n_nodes, ar_model_order, ...
	n_samples_cross_coup, n_samples_noise_corr);

% the same for all measures
blank_struct.model		= model_specific_filename;
blank_struct.n_nodes		= n_nodes;
blank_struct.spectral_radius	= spectral_radius;
blank_struct.cross_coup_range	= cross_coup_range;
blank_struct.rmi_range		= rmi_range;
blank_struct.coupling_mag	= cnzeros;
blank_struct.noise_corr		= nncnzeros;
blank_struct.coupling_matrices	= nnmcnzeros;

% measure-dependent
switch measure

case 'multi_info'
	blank_struct.MultiInfo			= cnzeros;

case 'dd_ce_co_info'
       	blank_struct.m_dim_ce_dd_co_info	= m_dim_ce_dd_co_info;

	if any(strcmp(dim_reduction,'pca'))
		blank_struct.DD_PCA		= cnzeros;
		blank_struct.ShannonCE_PCA	= cnzeros;
		blank_struct.CoInfoPCA		= cnzeros;
	end

	if any(strcmp(dim_reduction,'grassmanian'))
		blank_struct.DDGrassMin		= cnzeros;
		blank_struct.DDGrassMean	= cnzeros;
		blank_struct.ShannonCEGrassMax	= cnzeros;
		blank_struct.ShannonCEGrassMean	= cnzeros;
		blank_struct.CoInfoGrassMin	= cnzeros;
		blank_struct.CoInfoGrassMax	= cnzeros;
		blank_struct.CoInfoGrassMean	= cnzeros;
	end

case 'phiid_measures_mmi'
	blank_struct.PhiID_DoubleRed_MMI	= cnzeros;
	blank_struct.PhiID_DoubleSyn_MMI	= cnzeros;
	blank_struct.PhiID_CE_MMI		= cnzeros;
	blank_struct.PhiID_DC_MMI		= cnzeros;
	blank_struct.PhiID_CD_MMI		= cnzeros;
	blank_struct.PhiID_UC_MMI		= cnzeros;
	blank_struct.PhiID_Syn_MMI		= cnzeros;
	blank_struct.PhiID_Transfer_MMI		= cnzeros;

case 'phiid_measures_ccs'
	blank_struct.PhiID_DoubleRed_CCS	= cnzeros;
	blank_struct.PhiID_DoubleSyn_CCS	= cnzeros;
	blank_struct.PhiID_CE_CCS		= cnzeros;
	blank_struct.PhiID_DC_CCS		= cnzeros;
	blank_struct.PhiID_CD_CCS		= cnzeros;
	blank_struct.PhiID_UC_CCS		= cnzeros;
	blank_struct.PhiID_Syn_CCS		= cnzeros;
	blank_struct.PhiID_Transfer_CCS	= cnzeros;

case 'integrated_info_measures'
	if any(strcmp(list_integrated_info_measures,'AverageCorrelation'   )), blank_struct.AverageCorr			= cnzeros; end
	if any(strcmp(list_integrated_info_measures,'IntegratedInformation')), blank_struct.IntegratedInfo		= cnzeros; end
	if any(strcmp(list_integrated_info_measures,'IntegratedInteraction')), blank_struct.IntegratedInteraction	= cnzeros; end
	if any(strcmp(list_integrated_info_measures,'DecoderIntegration'   )), blank_struct.DecoderIntegration		= cnzeros; end
	if any(strcmp(list_integrated_info_measures,'CausalDensity'        )), blank_struct.CausalDensity		= cnzeros; end
	if any(strcmp(list_integrated_info_measures,'IntegratedSynergy'    )), blank_struct.IntegratedSynergy		= cnzeros; end
	if any(strcmp(list_integrated_info_measures,'TimeDelayedMutualInfo')), blank_struct.TimeDelayedMI		= cnzeros; end

otherwise error('Unknown measure ''%s''',measure);
end % switch measure
	
results = repmat(blank_struct, n_cross_coup, n_rmi); % preallocate

w = whos('results');
fprintf('Pre-allocated size of ''results'' = %.4f MB\n',w.bytes/1000/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end memory pre-allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logical indices of on/off-diagonal matrix elements for VAR coefficients sequence
arondiag  = repmat(diag(true(n_nodes,1)),1,1,ar_model_order);
aroffdiag = ~arondiag;

all_cross_coup_mag = 2.^linspace(cross_coup_range(1),cross_coup_range(2), n_cross_coup); % norm values (exponential scale)
all_rmi   = 2.^linspace(rmi_range(1), rmi_range(2),  n_rmi  ); % rmi values  (exponential scale)

for i = 1:n_cross_coup;
	fprintf('\nNORM %d of %d = %g\n', i, n_cross_coup, all_cross_coup_mag(i));

	for j = 1:n_rmi;
		fprintf('\n\tRMI %d of %d = %g\n', j, n_rmi, all_rmi(j));
		fprintf('\nNORM %d of %d = %g\n', i, n_cross_coup, all_cross_coup_mag(i));
		st = tic;

		for k = 1:n_samples_cross_coup;
			fprintf('\t\tcouplings sample %d of %d\n', k, n_samples_cross_coup);

			for p = 1:n_samples_noise_corr;
				fprintf('\t\t\tnoise corr sample %d of %d\n', p, n_samples_noise_corr);

				%------------------------------------------------------------------------
				% COUPLING & NOISE CORRELATION MATRICES
				%------------------------------------------------------------------------

				% TODO - Matrix 1.0e+07 norm does not seem to me to define "coupling strength"
				% in any meaningful sense. Furthermore, scaling a VAR coefficients sequence
				% actually fails for 1-lag VAR matrices, since the scaling is simply reversed
				% by enforcing the spectral radius!
				%
				% So here we construct a random VAR coefficients sequence (with spectral norm 1
				% to set a "baseline"); we then scale just the *off-diagonal* elements before
				% imposing the specified spectral radius.

				% construct VAR coefficients array of specified lags with spectral norm = 1
				% var_rand() does set a specific seed, no other seeds should therefore be set,
				% as that leads to repetitions of starting points and hence duplicate results
				var_coeff_seq = var_rand(n_nodes, ar_model_order, 1);

				% scale all off-diagonal elements
				var_coeff_seq(aroffdiag) = all_cross_coup_mag(i)*var_coeff_seq(aroffdiag);

				% enforce spectral norm as specified
				var_coeff_seq = specnorm(var_coeff_seq, spectral_radius);

				% TODO - It's still not clear to me what "coupling strength" is really
				% supposed to mean. In the *functional* sense, multi-information could in
				% fact be interpreted as a coupling strength, or for directed coupling,
				% "global" GC. It seems that what we want is some "network-structural"
				% notion of coupling strength, but what might that be? Below is a plausible
				% implementation, which identifies strong coupling if the off-diagonal
				% entries of the VAR coefficient matrices are on average much larger than
				% the the on-diagonal entries.
				%
				% Specifically, we define coupling strength as the ratio of the root-mean
				% -square of off-diagonal entries to the root-mean-square of all entries;
				% this will be zero if off-diagonal entries are zero (no coupling), and
				% approach 1 for off-diag entries >> on-diag entries (strong coupling).
				% This is pretty crude, though, and the variance will be high for given
				% spectral radius and scaling.

				results(i,j).coupling_mag(k,p) = sqrt(mean(var_coeff_seq(aroffdiag).^2)/mean(var_coeff_seq(:).^2));

				% generate random noise correlation matrices
				% note: -rmi sets rmi relative to uniform correlation matrix
				% at the given diemnsion
				noise_corr = corr_rand(n_nodes, -all_rmi(j));

				results(i,j).noise_corr(:,:,k,p) = noise_corr;

				results(i,j).coupling_matrices(:,:,:,k,p) = var_coeff_seq;

				% calculate autocovariance sequence (up to 1 lag only)
				autocov_seq = var_to_autocov(var_coeff_seq, noise_corr, -1); % -1 means *exactly* one lag!

				covmat  = autocov_seq(:,:,1);                % unlagged covariance
				ac1mat  = autocov_seq(:,:,2);                % 1-lag autocovariance
				ac01mat = [covmat, ac1mat; ac1mat', covmat]; % 0,1-lag autocovariance (block-Toeplitz matrix)

				switch measure

				case 'multi_info'

					%------------------------------------------------------------------------
					% MULTI-INFORMATION
					%------------------------------------------------------------------------

					results(i,j).MultiInfo(k,p) = sum(log(diag(covmat))) - logdet(covmat); % process multi-information (nats)

					fprintf('\n\tMultiInfo = %g\n', results(i,j).MultiInfo(k,p));

				case 'dd_ce_co_info'

					%------------------------------------------------------------------------
					% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
					%------------------------------------------------------------------------

					% convert from VAR to ISS
					[A,C,K] = var_to_ss(var_coeff_seq, noise_corr);

					% pre-compute state prediction error covariance matrices
					[~,P] = iss2ce_precomp(A, C, K, noise_corr);

					% calculate DD & CE
					for z = 1:length(dim_reduction)

						if strcmp(dim_reduction{z}, 'pca')

							[~, L]				= pcacov(covmat);
							[CI, DD]			= iss2ce(L, A, C, K, noise_corr, covmat, P);
							CE				= CI - DD;
							results(i,j).DD_PCA(k,p)	= DD;
							results(i,j).ShannonCE_PCA(k,p)	= CE;
							results(i,j).CoInfoPCA(k,p)	= CI;

							fprintf('\n\tCoInfoPCA = %g\n', results(i,j).CoInfoPCA(k,p));

						elseif strcmp(dim_reduction{z}, 'grassmanian')

							L	= rand_orthonormal(n_nodes, m_dim_ce_dd_co_info, n_samples_dd_ce_co_info);
							CI	= zeros(n_samples_dd_ce_co_info,1);
							DD	= zeros(n_samples_dd_ce_co_info,1);

							for s = 1:n_samples_dd_ce_co_info
								[CI(s),DD(s)] = iss2ce(L(:,:,s), A, C, K, noise_corr, covmat, P);
							end

							CE					= CI - DD;
							results(i,j).DDGrassMin(k,p)		= min(DD);
							results(i,j).DDGrassMean(k,p)		= mean(DD);
							results(i,j).ShannonCEGrassMax(k,p)	= max(CE);
							results(i,j).ShannonCEGrassMean(k,p)	= mean(CE);
							results(i,j).CoInfoGrassMin(k,p)	= min(CI);
							results(i,j).CoInfoGrassMax(k,p)	= max(CI);
							results(i,j).CoInfoGrassMean(k,p)	= mean(CI);

							fprintf('\n\tCoInfoGrassMax = %g\n', results(i,j).CoInfoGrassMax(k,p));

						end

					end

				case 'phiid_measures_mmi'

                                        %------------------------------------------------------------------------
                                        % PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, 
						    % CAUSAL DECOUPLING USING MMI
                                        %------------------------------------------------------------------------

                                        PhiID_Atoms_MMI		= PhiIDFullAnalytical(ac01mat, 'MMI');

                                        PhiID_CE_MMI		= PhiID_Atoms_MMI.str + PhiID_Atoms_MMI.stx + ...
							    PhiID_Atoms_MMI.sty + PhiID_Atoms_MMI.sts;
                                        PhiID_DC_MMI		= PhiID_Atoms_MMI.str + PhiID_Atoms_MMI.stx + ...
							    PhiID_Atoms_MMI.sty;
                                        PhiID_CD_MMI		= PhiID_CE_MMI - PhiID_DC_MMI;
                                        PhiID_DoubleRed_MMI	= PhiID_Atoms_MMI.rtr;
                                        PhiID_DoubleSyn_MMI	= PhiID_Atoms_MMI.sts;
                                        PhiID_UC_MMI		= PhiID_Atoms_MMI.rts + PhiID_Atoms_MMI.xts + ...
							    PhiID_Atoms_MMI.yts;
                                        PhiID_Syn_MMI		= PhiID_UC_MMI + PhiID_DC_MMI + PhiID_DoubleSyn_MMI;
                                        PhiID_Transfer_MMI	= PhiID_Atoms_MMI.xty + PhiID_Atoms_MMI.ytx;

                                        %results_phiid_atoms(i,j).PhiID_Atoms_MMI(k,p) = PhiID_Atoms_MMI;
                                        results(i,j).PhiID_CE_MMI(k,p)		= PhiID_CE_MMI;
                                        results(i,j).PhiID_DC_MMI(k,p)		= PhiID_DC_MMI;
						    results(i,j).PhiID_CD_MMI(k,p)		= PhiID_CD_MMI;
                                        results(i,j).PhiID_DoubleRed_MMI(k,p)	= PhiID_DoubleRed_MMI;
                                        results(i,j).PhiID_DoubleSyn_MMI(k,p)	= PhiID_DoubleSyn_MMI;
                                        results(i,j).PhiID_UC_MMI(k,p)		= PhiID_UC_MMI;
                                        results(i,j).PhiID_Syn_MMI(k,p)		= PhiID_Syn_MMI;
                                        results(i,j).PhiID_Transfer_MMI(k,p)	= PhiID_Transfer_MMI;

                                        fprintf('\n\tPhiID_CE_MMI = %g\n', PhiID_CE_MMI);
						    
				case 'phiid_measures_ccs'

                                        %------------------------------------------------------------------------
                                        % PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, 
						    % CAUSAL DECOUPLING USING CCS
                                        %------------------------------------------------------------------------

						    [micro corr_micro] = sim_mvar_model_with_lagged_cov(n_datapoints_phiid_ccs, ...
							    ac01mat, var_coeff_seq, ar_model_order);
						    
						    PhiID_Atoms_CCS	= PhiIDFullContinuousData(micro, ar_model_order, 'CCS', ...
							    'Gaussian');
						    %PhiID_Atoms_CCS	= PhiIDFull_Analytical(ac01mat, 'CCS');
	
						    PhiID_CE_CCS	= PhiID_Atoms_CCS.str + PhiID_Atoms_CCS.stx + ...
							    PhiID_Atoms_CCS.sty + PhiID_Atoms_CCS.sts;
						    PhiID_DC_CCS	= PhiID_Atoms_CCS.str + PhiID_Atoms_CCS.stx + ...
							    PhiID_Atoms_CCS.sty;
						    PhiID_CD_CCS	= PhiID_CE_CCS - PhiID_DC_CCS;
						    PhiID_DoubleRed_CCS	= PhiID_Atoms_CCS.rtr;
						    PhiID_DoubleSyn_CCS	= PhiID_Atoms_CCS.sts;
						    PhiID_UC_CCS	= PhiID_Atoms_CCS.rts + PhiID_Atoms_CCS.xts + ...
							    PhiID_Atoms_CCS.yts;
						    PhiID_Syn_CCS	= PhiID_UC_CCS + PhiID_DC_CCS + PhiID_DoubleSyn_CCS;
						    PhiID_Transfer_CCS	= PhiID_Atoms_CCS.xty + PhiID_Atoms_CCS.ytx;
	
						    %results_phiid_atoms(i,j).PhiID_Atoms_CCS(k,p) = PhiID_Atoms_CCS;
						    results(i,j).PhiID_CE_CCS(k,p)		= PhiID_CE_CCS;
						    results(i,j).PhiID_DC_CCS(k,p)		= PhiID_DC_CCS;
						    results(i,j).PhiID_CD_CCS(k,p)		= PhiID_CD_CCS;
						    results(i,j).PhiID_DoubleRed_CCS(k,p)	= PhiID_DoubleRed_CCS;
						    results(i,j).PhiID_DoubleSyn_CCS(k,p)	= PhiID_DoubleSyn_CCS;
						    results(i,j).PhiID_UC_CCS(k,p)		= PhiID_UC_CCS;
						    results(i,j).PhiID_Syn_CCS(k,p)		= PhiID_Syn_CCS;
						    results(i,j).PhiID_Transfer_CCS(k,p)	= PhiID_Transfer_CCS;
						    
						    fprintf('\n\tPhiID_CE_CCS = %g\n', PhiID_CE_CCS);

				case 'integrated_info_measures'

					%------------------------------------------------------------------------
					% INTEGRATED INFORMATION MEASURES
					%------------------------------------------------------------------------

					% makeLaggedCovariance calculates the covariance, but with ones in the diagonals,
					% and values ranging from -1 to 1 in the off-diagonals; this full time-lagged
					% covariance matrix will be different from the one above, as the one
					% above normalizes coupling matrices with a decay factor using specnorm(),
					% therefore, we'll use the one above to keep things consistent

					for z = 1:length(list_integrated_info_measures);

						% name template to instantiate JIDT calculators
						class_template = 'infodynamics.measures.continuous.gaussian.%sCalculatorGaussian';

						% the following uses only the full time-lagged covariance matrix to compute the measure
						if strcmp(list_integrated_info_measures{z}, 'AverageCorrelation')
							diag_full_time_lagged_cov = diag(diag(ac01mat));

							slS	= (diag_full_time_lagged_cov^-0.5) * ...
								ac01mat * (diag_full_time_lagged_cov^-0.5);

							measure_value = (sum(abs(slS(:))) - ...
								sum(diag(abs(slS))))/(n_nodes^2-n_nodes);
						else
							% TODO: My understanding is that the following code wants a *single-lag* "coupling matrix".
							% The next line is (probably?) okay if a single lag is specified (ar_model_order = 1), not
							% so sure about ar_model_order > 1 ... another problem in that case is that 'ac01mat' will
							% no longer be correct as calculated previously.
							%
							% TODO next line commented out because not used here?
							% coupling_matrix = var_coeff_seq(:,:,1); % set to lag-1 coefficients matrix
							try
								calc = javaObject(sprintf(class_template, list_integrated_info_measures{z}));
								calc.initialise(n_nodes); % TODO - is this number of nodes???
								calc.setLaggedCovariance(ac01mat);
								measure_value = calc.computeAverageLocalOfObservations();
							catch
								measure_value = NaN;
							end
						end

						switch list_integrated_info_measures{z}
							case 'AverageCorrelation',    results(i,j).AverageCorr(k,p)		= measure_value;
							case 'IntegratedInformation', results(i,j).IntegratedInfo(k,p)		= measure_value;
							case 'IntegratedInteraction', results(i,j).IntegratedInteraction(k,p)	= measure_value;
							case 'DecoderIntegration',    results(i,j).DecoderIntegration(k,p)	= measure_value;
							case 'CausalDensity',         results(i,j).CausalDensity(k,p)		= measure_value;
							case 'IntegratedSynergy',     results(i,j).IntegratedSynergy(k,p)	= measure_value;
							case 'TimeDelayedMutualInfo', results(i,j).TimeDelayedMI(k,p)       	= measure_value;
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
		fprintf('\t\tcoupling/noise loop completed in %g seconds; size of ''results'' = %.4f MB\n', et, w.bytes/1000/1000);

	end

end

% save results

pathout_data_measures  = fullfile(basedir, 'results', 'analyses', 'multivar_autoreg', 'rand_cross_coup_rmi', model_specific_path);

% row & column values for table
all_cross_coup_mag_str = {};
for t = 1:length(all_cross_coup_mag)
        all_cross_coup_mag_str{t} = num2str(all_cross_coup_mag(t));
end

all_rmi_str = {};
for e = 1:length(all_rmi)
        all_rmi_str{e} = num2str(all_rmi(e));
end


if strcmp(measure, 'phiid_based_measures_mmi') || strcmp(measure, 'phiid_based_measures_ccs')

	if strcmp(measure, 'phiid_based_measures_mmi')
		phiid_fieldnames = {'PhiID_DoubleRed_MMI', 'PhiID_DoubleSyn_MMI', ...
 		'PhiID_CE_MMI', 'PhiID_DC_MMI', 'PhiID_CD_MMI', 'PhiID_UC_MMI', ...
 		'PhiID_Syn_MMI', 'PhiID_Transfer_MMI'};

 	elseif strcmp(measure, 'phiid_based_measures_ccs')
 		phiid_fieldnames = {'PhiID_DoubleRed_CCS', 'PhiID_DoubleSyn_CCS', ...
 		'PhiID_CE_CCS', 'PhiID_DC_CCS', 'PhiID_CD_CCS', 'PhiID_UC_CCS', ...
 		'PhiID_Syn_CCS', 'PhiID_Transfer_CCS'};

 	end
 	
 	for h = 1:length(phiid_fieldnames)
		for i=1:size(results,1)
			for j=1:size(results,2)
			
			results_phiid(i,j).n_nodes	         = results(i,j).n_nodes;
			results_phiid(i,j).model		 = results(i,j).model;
			results_phiid(i,j).n_nodes		 = results(i,j).n_nodes;
			results_phiid(i,j).spectral_radius	 = results(i,j).spectral_radius;
			results_phiid(i,j).cross_coup_range	 = results(i,j).cross_coup_range;
			results_phiid(i,j).rmi_range		 = results(i,j).rmi_range;
			results_phiid(i,j).coupling_mag		 = results(i,j).coupling_mag;
			results_phiid(i,j).noise_corr	         = results(i,j).noise_corr;
			results_phiid(i,j).coupling_matrices	 = results(i,j).coupling_matrices;
			results_phiid(i,j).(phiid_fieldnames{h}) = results(i,j).(phiid_fieldnames{h});

			end
		end
	
		% save as variable
		results_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' results_filename '_' phiid_fieldnames{h} '.mat']);
		fprintf('\nSaving results file ''%s'' ...', results_filepath);
		save(results_filepath,'results', '-v7.3');
		fprintf(' done\n');
		
		% as table
		results_table_filename = [model_specific_filename '_' results_filename '_' phiid_fieldnames{h} '_table'];
		results_table = array2table(results_phiid, 'RowNames', all_cross_coup_mag_str, 'VariableNames', all_rmi_str);
		results_table_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' results_table_filename '.mat']);
		fprintf('\nSaving results table file ''%s'' ...', results_table_filepath);
		save(results_table_filepath, 'results_table', '-v7.3');
		fprintf(' d1.0e+07 one\n');
		
	end

else 
		
	% save as variable

	results_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' results_filename '.mat']);
	fprintf('\nSaving results file ''%s'' ...', results_filepath);
	save(results_filepath,'results', '-v7.3');
	fprintf(' done\n');

	% save as table

	results_table_filename = [results_filename '_table'];
	results_table = array2table(results, 'RowNames', all_cross_coup_mag_str, 'VariableNames', all_rmi_str);

	results_table_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' results_table_filename '.mat']);
	fprintf('\nSaving results table file ''%s'' ...', results_table_filepath);
	save(results_table_filepath, 'results_table', '-v7.3');
	fprintf(' d1.0e+07 one\n');
end

