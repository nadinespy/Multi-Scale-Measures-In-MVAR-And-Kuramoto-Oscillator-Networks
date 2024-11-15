i%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% measures_random_couplings_lionel.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this script, set 'basedir' and 'measure', and run the appropriate parameter script
%
% E.g., to run interactively:
%
% >> basedir = '/its/home/ns508';
% >> measure = 'DD_CE_CO_INFO';
% >> mrc_parameters;
% >> measures_random_couplings_lionel;
%
% To run in batch mode, See Bash script 'bash_measures_random_couplings_lionel.sh'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display parameters

fprintf('\nParameters:\n-----------\n\n');
fprintf('current directory     : %s\n',      pwd);
fprintf('base directory        : %s\n',      basedir);
fprintf('measure               : %s\n',      measure);
fprintf('n_nodes               : %d\n',      n_nodes);
if exist('m_dim')
	fprintf('m_dim                 : %d\n',      m_dim);
end
if exist('dim_reduction')
	fprintf('dim_reduction         : %s\n',      strjoin(dim_reduction,', '));
end
fprintf('ar_model_order        : %d\n',      ar_model_order);
fprintf('spectral_radius       : %g\n',      spectral_radius);
fprintf('n_samples_couplings   : %d\n',      n_samples_couplings);
fprintf('n_samples_noise_corrs : %d\n',      n_samples_noise_corrs);
fprintf('n_samples_ddce        : %d\n',      n_samples_noise_corrs);
fprintf('norm_range            : %g - %g\n', norm_range(1),norm_range(2));
fprintf('n_norms               : %d\n',      n_norms);
fprintf('rmi_range             : %g - %g\n', rmi_range(1),rmi_range(2));
fprintf('n_rmi                 : %d\n',      n_rmi);
fprintf('results_filename      : %s\n',      results_filename);
fprintf('\n');

% run SSDI startup (also calls MVGC2 startup)
run(fullfile(basedir,'packages_and_code_repos','ssdi','startup'));

addpath(fullfile(basedir,'packages_and_code_repos','ReconcilingEmergences-master'));
addpath(fullfile(basedir,'code','common'));
addpath(fullfile(basedir,'code','analytical','functions'));
addpath(fullfile(basedir,'code','analytical','scripts'));
addpath(fullfile(basedir,'code','analytical','scripts','measures_random_couplings'));

% {
javaaddpath(fullfile(basedir,'/infodynamics.jar'));
javaaddpath(fullfile(basedir,'code','analytical','scripts','measures_random_couplings','infodynamics.jar'));
javaaddpath(fullfile(basedir,'code','analytical','scripts','measures_random_couplings','commons-math3-3.5.jar'));
javaaddpath(fullfile(basedir,'code','common','infodynamics.jar'));
javaaddpath(fullfile(basedir,'code','common','commons-math3-3.5.jar'));

try
	javaaddpath(fullfile(basedir,'code','analytical','scripts','measures_random_couplings','infodynamics.jar'));
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
cnzeros    = zeros(n_samples_couplings, n_samples_noise_corrs);
nncnzeros  = zeros(n_nodes, n_nodes, n_samples_couplings, n_samples_noise_corrs);
nnmcnzeros = zeros(n_nodes, n_nodes, ar_model_order, n_samples_couplings, n_samples_noise_corrs);

% the same for all measures
blank_struct.model	       = model_specific_filename;
blank_struct.n_nodes           = n_nodes;
blank_struct.spectral_radius   = spectral_radius;
blank_struct.cross_coup_range  = norm_range;
blank_struct.rmi_range         = rmi_range;
blank_struct.coupling_mag      = cnzeros;
blank_struct.noise_corr        = nncnzeros;
blank_struct.coupling_matrices = nnmcnzeros;

% measure-dependent
switch measure
case 'multi_info'
	blank_struct.multi_info = cnzeros;

case 'DD_CE_CO_INFO'

       	blank_struct.m_dim           = m_dim;

	if any(strcmp(dim_reduction,'pca'))
		blank_struct.DD_pca                = cnzeros;
		blank_struct.ShannonCE_pca         = cnzeros;
		blank_struct.co_info_pca           = cnzeros;
	end
	if any(strcmp(dim_reduction,'grassmanian'))
		blank_struct.DD_grassmanian        = cnzeros;
		blank_struct.ShannonCE_grassmanian = cnzeros;
		blank_struct.co_info_grassmanian   = cnzeros;
	end

case 'phiid_based_measures_mmi'
	blank_struct.phiidDoubleRed_MMI = cnzeros;
	blank_struct.phiidDoubleSyn_MMI = cnzeros;
	blank_struct.phiidCE_MMI        = cnzeros;
	blank_struct.phiidDC_MMI        = cnzeros;
	blank_struct.phiidCD_MMI        = cnzeros;
	blank_struct.phiidUC_MMI        = cnzeros;
	blank_struct.phiidSyn_MMI       = cnzeros;
	blank_struct.phiidTransfer_MMI  = cnzeros;

case 'phiid_based_measures_ccs'
	blank_struct.phiidDoubleRed_CCS = cnzeros;
	blank_struct.phiidDoubleSyn_CCS = cnzeros;
	blank_struct.phiidCE_CCS        = cnzeros;
	blank_struct.phiidDC_CCS        = cnzeros;
	blank_struct.phiidCD_CCS        = cnzeros;
	blank_struct.phiidUC_CCS        = cnzeros;
	blank_struct.phiidSyn_CCS       = cnzeros;
	blank_struct.phiidTransfer_CCS  = cnzeros;

case 'integrated_info_measures'
	if any(strcmp(calcNames,'AverageCorrelation'   )), blank_struct.average_corr           = cnzeros; end
	if any(strcmp(calcNames,'IntegratedInformation')), blank_struct.integrated_info        = cnzeros; end
	if any(strcmp(calcNames,'IntegratedInteraction')), blank_struct.integrated_interaction = cnzeros; end
	if any(strcmp(calcNames,'DecoderIntegration'   )), blank_struct.decoder_integration    = cnzeros; end
	if any(strcmp(calcNames,'CausalDensity'        )), blank_struct.causal_density         = cnzeros; end
	if any(strcmp(calcNames,'IntegratedSynergy'    )), blank_struct.integrated_synergy     = cnzeros; end
	if any(strcmp(calcNames,'TimeDelayedMutualInfo')), blank_struct.time_delayed_mi        = cnzeros; end

otherwise error('Unknown measure ''%s''',measure);
end % switch measure

results = repmat(blank_struct, n_norms, n_rmi); % preallocate

w = whos('results');
fprintf('Pre-allocated size of ''results'' = %.4f MB\n',w.bytes/1000/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end memory pre-allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logical indices of on/off-diagonal matrix elements for VAR coefficients sequence
arondiag  = repmat(diag(true(n_nodes,1)),1,1,ar_model_order);
aroffdiag = ~arondiag;

all_norms = 2.^linspace(norm_range(1),norm_range(2), n_norms); % norm values (exponential scale)
all_rmi   = 2.^linspace(rmi_range(1), rmi_range(2),  n_rmi  ); % rmi values  (exponential scale)

for i = 1:n_norms;
	fprintf('\nNORM %d of %d = %g\n',i,n_norms,all_norms(i));

	for j = 1:n_rmi;
		fprintf('\n\tRMI %d of %d = %g\n',j,n_rmi,all_rmi(j));
		fprintf('\nNORM %d of %d = %g\n',i,n_norms,all_norms(i));
		st = tic;

		for k = 1:n_samples_couplings;
			fprintf('\t\tcouplings sample %d of %d\n',k,n_samples_couplings);

			for p = 1:n_samples_noise_corrs;
				fprintf('\t\t\tnoise corr sample %d of %d\n',p,n_samples_noise_corrs);

				%------------------------------------------------------------------------
				% COUPLING & NOISE CORRELATION MATRICES
				%------------------------------------------------------------------------

				% TODO - Matrix norm does not seem to me to define "coupling strength"
				% in any meaningful sense. Furthermore, scaling a VAR coefficients sequence
				% actually fails for 1-lag VAR matrices, since the scaling is simply reversed
				% by enforcing the spectral radius!
				%
				% So here we construct a random VAR coefficients sequence (with spectral norm 1
				% to set a "baseline"); we then scale just the *off-diagonal* elements before
				% imposing the specified spectral radius.

				% construct VAR coefficients array of specified lags with spectral norm = 1
				var_coeff_seq = var_rand(n_nodes,ar_model_order,1);

				% scale all off-diagonal elements
				var_coeff_seq(aroffdiag) = all_norms(i)*var_coeff_seq(aroffdiag);

				% enforce spectral norm as specified
				var_coeff_seq = specnorm(var_coeff_seq,spectral_radius);

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
				% at the given dimension
				noise_corr = corr_rand(n_nodes,-all_rmi(j));

				results(i,j).noise_corr(:,:,k,p) = noise_corr;

				results(i,j).coupling_matrices(:,:,:,k,p) = var_coeff_seq;

				% calculate autocovariance sequence (up to 1 lag only)
				autocov_seq = var_to_autocov(var_coeff_seq,noise_corr,-1); % -1 means *exactly* one lag!

				covmat  = autocov_seq(:,:,1);                % unlagged covariance
				ac1mat  = autocov_seq(:,:,2);                % 1-lag autocovariance
				ac01mat = [covmat, ac1mat; ac1mat', covmat]; % 0,1-lag autocovariance (block-Toeplitz matrix)

				switch measure

				case 'multi_info'

					%------------------------------------------------------------------------
					% MULTI-INFORMATION
					%------------------------------------------------------------------------

					results(i,j).multi_info(k,p) = sum(log(diag(covmat))) - logdet(covmat); % process multi-information (nats)
					
					fprintf('\n\tmulti_info = %g\n', results(i,j).multi_info(k,p));

				case 'DD_CE_CO_INFO'

					%------------------------------------------------------------------------
					% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
					%------------------------------------------------------------------------

					% convert from VAR to ISS
					[A,C,K] = var_to_ss(var_coeff_seq,noise_corr);

					% pre-compute state prediction error covariance matrices
					[~,P] = iss2ce_precomp(A,C,K,noise_corr);

					% calculate DD & CE
					for z = 1:length(dim_reduction)

						if strcmp(dim_reduction{z}, 'pca')

							[~,L] = pcacov(covmat);
							[CI,DD] = iss2ce(L,A,C,K,noise_corr,covmat,P);
							CE = CI - DD;
							results(i,j).DD_pca(k,p)        = DD;
							results(i,j).ShannonCE_pca(k,p) = CE;
							results(i,j).co_info_pca(k,p)   = CI;

						elseif strcmp(dim_reduction{z}, 'grassmanian')

							L = rand_orthonormal(n_nodes,m_dim,n_samples_ddce);
							CI = zeros(n_samples_ddce,1);
							DD = zeros(n_samples_ddce,1);
							for s = 1:n_samples_ddce
								[CI(s),DD(s)] = iss2ce(L(:,:,s),A,C,K,noise_corr,covmat,P);
							end
							CE = CI - DD;
							results(i,j).DD_grassmanian(k,p)        = min(DD);
							results(i,j).ShannonCE_grassmanian(k,p) = max(CE);
							results(i,j).co_info_grassmanian(k,p)   = max(CI);

						end

					end

				case 'phiid_based_measures_mmi'

					%------------------------------------------------------------------------
					% PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, CAUSAL DECOUPLING
					%------------------------------------------------------------------------

					% MMI
					phiid_atoms_CCE = PhiIDFull_Analytical(ac01mat, 'MMI');

					phiidCE_MMI        = phiid_atoms.str + phiid_atoms.stx + phiid_atoms.sty + phiid_atoms.sts;
					phiidDC_MMI        = phiid_atoms.str + phiid_atoms.stx + phiid_atoms.sty;
					phiidCD_MMI        = phiidCE - phiidDC;
					phiidDoubleRed_MMI = phiid_atoms.rtr;
					phiidDoubleSyn_MMI = phiid_atoms.sts;
					phiidUC_MMI        = phiid_atoms.rts + phiid_atoms.xts + phiid_atoms.yts;
					phiidSyn_MMI       = phiidUC + phiidDC + phiidDoubleSyn;
					phiidTransfer_MMI  = phiid_atoms.xty + phiid_atoms.ytx;

					results(i,j).phiid_atoms_MMI(k,p)    = phiid_atoms_MMI;
					results(i,j).phiidCE_MMI(k,p)        = phiidCE_MMI;
					results(i,j).phiidDC_MMI(k,p)        = phiidDC_MMI;
					results(i,j).phiidCD_MMI(k,p)        = phiidCD_MMI;
					results(i,j).phiidDoubleRed_MMI(k,p) = phiidDoubleRed_MMI;
					results(i,j).phiidDoubleSyn_MMI(k,p) = phiidDoubleSyn_MMI;
					results(i,j).phiidUC_MMI(k,p)        = phiidUC_MMI;
					results(i,j).phiidSyn_MMI(k,p)       = phiidSyn_MMI;
					results(i,j).phiidTransfer_MMI(k,p)  = phiidTransfer_MMI;
					
					fprintf('\n\tphiidCE = %g\n', phiidCE_MMI);

				case 'phiid_based_measures_ccs'

					%------------------------------------------------------------------------
					% PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, CAUSAL DECOUPLING
					%------------------------------------------------------------------------
					
					% simulate time-series using full time-lagged covariance matrix
					time_length = 2000;
					[X corr_X] = sim_mvar_network_with_lagged_cov(time_length, ac01mat, ar_model_order)

					% CCS
					phiid_atoms_CCE = PhiIDFull_Analytical(ac01mat, 'CCS');

					phiidCE_CCS        = phiid_atoms.str + phiid_atoms.stx + phiid_atoms.sty + phiid_atoms.sts;
					phiidDC_CCS        = phiid_atoms.str + phiid_atoms.stx + phiid_atoms.sty;
					phiidCD_CCS        = phiidCE - phiidDC;
					phiidDoubleRed_CCS = phiid_atoms.rtr;
					phiidDoubleSyn_CCS = phiid_atoms.sts;
					phiidUC_CCS        = phiid_atoms.rts + phiid_atoms.xts + phiid_atoms.yts;
					phiidSyn_CSS       = phiidUC + phiidDC + phiidDoubleSyn;
					phiidTransfer_CCS  = phiid_atoms.xty + phiid_atoms.ytx;

					results(i,j).phiid_atoms_CCS(k,p)    = phiid_atoms_CCS;
					results(i,j).phiidCE_CCS(k,p)        = phiidCE_CCS;
					results(i,j).phiidDC_CCS(k,p)        = phiidDC_CCS;
					results(i,j).phiidCD_CCS(k,p)        = phiidCD_CCS;
					results(i,j).phiidDoubleRed_CCS(k,p) = phiidDoubleRed_CCS;
					results(i,j).phiidDoubleSyn_CCS(k,p) = phiidDoubleSyn_CCS;
					results(i,j).phiidUC_CCS(k,p)        = phiidUC_CCS;
					results(i,j).phiidSyn_CCS(k,p)       = phiidSyn_CCS;
					results(i,j).phiidTransfer_CCS(k,p)  = phiidTransfer_CCS;

					fprintf('\n\tphiidCE_CCS = %g\n', phiidCE_CCS);

				case 'integrated_info_measures'

					%------------------------------------------------------------------------
					% INTEGRATED INFORMATION MEASURES
					%------------------------------------------------------------------------

					% makeLaggedCovariance calculates the covariance, but with ones in the diagonals,
					% and values ranging from -1 to 1 in the off-diagonals; this full time-lagged
					% covariance matrix will be different from the one above, as the one
					% above normalizes coupling matrices with a decay factor using specnorm(),
					% therefore, we'll use the one above to keep things consistent

					for z = 1:length(calcNames);

						% name template to instantiate JIDT calculators
						class_template = 'infodynamics.measures.continuous.gaussian.%sCalculatorGaussian';

						% the following uses only the full time-lagged covariance matrix to compute the measure
						if strcmp(calcNames{z}, 'AverageCorrelation')
							diag_full_time_lagged_cov = diag(diag(ac01mat));
							slS = (diag_full_time_lagged_cov^-0.5) * ac01mat * (diag_full_time_lagged_cov^-0.5);
							measure_value = (sum(abs(slS(:))) - sum(diag(abs(slS))))/(n_nodes^2-n_nodes);
						else
							% TODO: My understanding is that the following code wants a *single-lag* "coupling matrix".
							% The next line is (probably?) okay if a single lag is specified (ar_model_order = 1), not
							% so sure about ar_model_order > 1 ... another problem in that case is that 'ac01mat' will
							% no longer be correct as calculated previously.
							%
							% TODO next line commented out because not used here?
							% coupling_matrix = var_coeff_seq(:,:,1); % set to lag-1 coefficients matrix
							try
								calc = javaObject(sprintf(class_template, calcNames{z}));
								calc.initialise(n_nodes); % TODO - is this number of nodes???
								calc.setLaggedCovariance(ac01mat);
								measure_value = calc.computeAverageLocalOfObservations();
							catch
								measure_value = NaN;
							end
						end

						switch calcNames{z}
							case 'AverageCorrelation',    results(i,j).average_corr(k,p)           = measure_value;
							case 'IntegratedInformation', results(i,j).integrated_info(k,p)        = measure_value;
							case 'IntegratedInteraction', results(i,j).integrated_interaction(k,p) = measure_value;
							case 'DecoderIntegration',    results(i,j).decoder_integration(k,p)    = measure_value;
							case 'CausalDensity',         results(i,j).causal_density(k,p)         = measure_value;
							case 'IntegratedSynergy',     results(i,j).integrated_synergy(k,p)     = measure_value;
							case 'TimeDelayedMutualInfo', results(i,j).time_delayed_mi(k,p)        = measure_value;
							otherwise, error('Unknown ''calcName''');
						end % switch calcNames{z}
						fprintf('\n\tmeasure value = %g\n', measure_value);
					end
				end % switch(measure
			end
		end

		w = whos('results');
		et = toc(st);
		fprintf('\t\tcoupling/noise loop completed in %g seconds; size of ''results'' = %.4f MB\n',et,w.bytes/1000/1000);

	end

end

% save results

pathout_data_measures  = fullfile(basedir,'results','analyses', model_specific_path,'analytical','measures_random_couplings');

% as variable

results_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' results_filename '.mat']);
fprintf('\nSaving results file ''%s'' ...',results_filepath);
save(results_filepath,'results', '-v7.3');
fprintf(' done\n');

% as table

all_norms_str = {};
for t = 1:length(all_norms)
	all_norms_str{t} = num2str(all_norms(t));
end

all_rmi_str = {};
for e = 1:length(all_rmi)
	all_rmi_str{e} = num2str(all_rmi(e));
end

restable_filename = [results_filename '_table'];
results_table = array2table(results, 'RowNames', all_norms_str, 'VariableNames', all_rmi_str);

restable_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' restable_filename '.mat']);
fprintf('\nSaving results table file ''%s'' ...',restable_filepath);
save(restable_filepath,'results_table', '-v7.3');
fprintf(' done\n');
