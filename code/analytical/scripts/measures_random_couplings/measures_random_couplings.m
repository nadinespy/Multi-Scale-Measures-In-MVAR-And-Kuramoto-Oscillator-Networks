clear all;
close all;
clear java;
clc;

cd '/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/analytical/scripts'
javaaddpath('/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations/code/analytical/scripts/infodynamics.jar');
directories = @get_mvar_directories;

n_nodes			= 8;							% number of variables in network
m_dim				= 1;							% dimension of macro variable
dim_reduction		= {'pca', 'grassmanian'};			% dim reduction method: can be 'pca' or 'grassmanian' 
time_lag_for_model	= 1;							% number of time-lags
n_samples_couplings	= 50;							% number of random samples for connectivity matrices
n_samples_noise_corrs	= 25;							% number of random samples for noise correlation matrices
seed				= 1;

network			= [num2str(n_nodes) 'mvar' ...		% model name
				  '_lag' num2str(time_lag_for_model)];
			  
% full set of possible emergence measures
% measures			= {'DD_CE_CO_INFO', ...				% measures are clustered: if choosing 'DE_CE', co-info, DD & CE are 
% 				   'phiid_based_measures_mmi', ...		% calculated; if choosing 'integrated_info_measures', all measures
% 				   'phiid_based_measures_ccs', ...        % from [calcNames] are calculated; if choosing
% 				   'integrated_info_measures', ...		% 'phiid_based_measures_method', a bunch of PhiID-based measures 
% 				   'multi_info'};					% given the method are calculated
							
measures = {'DD_CE_CO_INFO', 'multi_info', 'integrated_info_measures'};:

% full set of possible integrated information measures to be calculated (alongside multi-info, DD, & CE)
% calcNames = {'IntegratedInformation', 'IntegratedInteraction', 'AverageCorrelation', ...
% 			'DecoderIntegration', 'CausalDensity', 'IntegratedSynergy', ...
% 			'TimeDelayedMutualInfo'};

calcNames = {'IntegratedInformation', 'AverageCorrelation', 'TimeDelayedMutualInfo'};
		
% idiosyncratic filename for this particular simulation run
filename_table = {'results_DD_CE_table'};
filename = {'results_DD_CE'};

rng(seed) % seed random number generator

% variation of residuals' mutual information and matrix norms
all_rmi		= linspace(0.0,1,10);		% array of rmi values 
all_norms		= linspace(0.01,1,10);		% array of norm-2 values
spectral_radius	= 0.9;

directories();
%% loop over coupling matrices and RMI values (including random sampling of connectivities & noise correlations)

set(0,'DefaultFigureWindowStyle','docked')

for i = 1:length(all_norms);
	
	for j = 1:length(all_rmi);
		
		for k = 1:n_samples_couplings;
			
			for p = 1:n_samples_noise_corrs;
				
				% define custom error message
				errID = 'DLYAP:NotDefinedOrNotUnique';
				msgtext = 'No or no unique solution to Lyapunov equation.';
				ME = MException(errID,msgtext);
				
				for l = 1:length(measures);
					
					% try to generate noise corrs/coupling matrices that give
					% unique DLYAP solutions, as long as the previous one was neither
					% unique nor existent
					count = 0;
					err_count = 0;
					while count == err_count
						
						try
							%------------------------------------------------------------------------
							% COUPLING & NOISE CORRELATION MATRICES
							%------------------------------------------------------------------------
							% create random 2-node matrices of the form [a b; c d] for norm value all_norms(i)
							coupling_matrix = randn(n_nodes);
							coupling_matrix = all_norms(i)*coupling_matrix/norm(coupling_matrix);
							
							% append to account for time-lag
							full_coupling_matrix = [];
							for g = 1:time_lag_for_model+1
								full_coupling_matrix = cat(3, full_coupling_matrix, coupling_matrix);
							end
							
							% normalize by spectral radius, and re-calculate spectral radius
							full_coupling_matrix = specnorm(full_coupling_matrix, spectral_radius);
							
							% spectral_radiuses(k,p) = max(abs(eig(full_coupling_matrix)));
							spectral_radiuses(k,p) = specnorm(full_coupling_matrix);
							
							% re-calculate the norm, as normalization by spectral radius will have changed it
							new_norms(k,p) = norm(full_coupling_matrix(:,:,1));
							
							% generate random noise correlation matrices
							noise_corr = corr_rand(n_nodes,all_rmi(j));
							
							noise_corrs{k,p}		= noise_corr;
							coupling_matrices{k,p}	= full_coupling_matrix;
							
							%------------------------------------------------------------------------
							% MULTI-INFORMATION
							%------------------------------------------------------------------------
							
							% process covariance matrix (solve DLYAP): cov_matrix(:,:,1) is the time-lagged covariance
							% matrix between t and t-1 etc, cov_matrix(:,:,2) the time-lagged covariance
							% matrix between t and t-2 etc.; cov_matrix(:,:,0) is the stationary covariance matrix (time-independent)
							
							cov_matrix		= var_to_autocov(full_coupling_matrix,noise_corr,time_lag_for_model);
							
							if strcmp(measures{l}, 'multi_info')
								
								total_cov(k,p)	= sum(log(diag(cov_matrix(:,:,1))));	% process total variance
								general_cov(k,p)	= log(det(cov_matrix(:,:,1)));		% process generalised variance (includes covariances)
								multi_info(k,p)	= total_cov(k,p) - general_cov(k,p);	% process multi-information (nats)
								
							elseif strcmp(measures{l}, 'DD_CE_CO_INFO')
								
								%------------------------------------------------------------------------
								% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
								%------------------------------------------------------------------------
								
								% run functions below without plotting and printing outputs (evalc() is a nasty solution)
								evalc('get_ssi_from_var(n_nodes, time_lag_for_model, coupling_matrix, full_coupling_matrix, spectral_radius, all_rmi(j), noise_corr)');
								
								% generate linear macro variables from the Grassmanian manifold
								% evalc('get_preoptimised_dd(mdim)');
								% evalc('get_optimised_dd(mdim)');
								
								% calculate DD & CE
								for z = 1:length(dim_reduction)
									
									evalc('[DD, CE] = get_ce_and_dd(m_dim, dim_reduction{z}, cov_matrix)');
									
									if strcmp(dim_reduction{z}, 'pca')
										
										all_DD_pca(k,p)		      = min(DD);
										all_ShannonCE_pca(k,p)	      = max(CE);
										all_co_info_pca(k,p)		= max(DD+CE);
										
									elseif strcmp(dim_reduction{z}, 'grassmanian')
										all_DD_grassmanian(k,p)		 = min(DD);
										all_ShannonCE_grassmanian(k,p) = max(CE);
										all_co_info_grassmanian(k,p)	 = min(DD) + max(CE);
										
									end
									
									clear DD;
									clear CE;
								end
								
							elseif strcmp(measures{l}, 'phiid_based_measures_mmi')
								
								%------------------------------------------------------------------------
								% PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, CAUSAL DECOUPLING
								%------------------------------------------------------------------------
								
								% construct full time-lagged covariance matrix
								
								% switch off-diagonal values in upper right quadrant
								if size(cov_matrix,1) < 3
									upper_right_quadrant = [cov_matrix(1,1,time_lag_for_model+1), cov_matrix(2,1,time_lag_for_model+1); ...
										cov_matrix(1,2,time_lag_for_model+1), cov_matrix(2,2,time_lag_for_model+1)];
								else
									upper_right_quadrant = cov_matrix(:,:,time_lag_for_model+1).';
								end
								
								% put all quadrants together
								full_time_lagged_cov = [cov_matrix(:,:,1), cov_matrix(:,:,time_lag_for_model+1); upper_right_quadrant, cov_matrix(:,:,1)];
								
								% MMI
								phiid_atoms_MMI		= PhiIDFull_Analytical(full_time_lagged_cov, 'MMI');
								
								phiidCE_MMI(k,p)		= phiid_atoms_MMI.str + phiid_atoms_MMI.stx + ...
									phiid_atoms_MMI.sty + phiid_atoms_MMI.sts;
								phiidDC_MMI(k,p)		= phiid_atoms_MMI.str + phiid_atoms_MMI.stx + ...
									phiid_atoms_MMI.sty;
								phiidCD_MMI(k,p)		= phiidCE_MMI(k,p) - phiidDC_MMI(k,p);
								phiidDoubleRed_MMI(k,p) = phiid_atoms_MMI.rtr;
								phiidDoubleSyn_MMI(k,p) = phiid_atoms_MMI.sts;
								
								phiidUC_MMI(k,p)		= phiid_atoms_MMI.rts + phiid_atoms_MMI.xts + phiid_atoms_MMI.yts;
								phiidSyn_MMI(k,p)		= phiidUC_MMI(k,p) + phiidDC_MMI(k,p) + phiidDoubleSyn_MMI(k,p);
								phiidTransfer_MMI(k,p)	= phiid_atoms_MMI.xty + phiid_atoms_MMI.ytx;
								
							elseif strcmp(measures{l}, 'phiid_based_measures_ccs')
								
								% construct full time-lagged covariance matrix
								
								% switch off-diagonal values in upper right quadrant
								if size(cov_matrix,1) < 3
									upper_right_quadrant = [cov_matrix(1,1,time_lag_for_model+1), cov_matrix(2,1,time_lag_for_model+1); ...
										cov_matrix(1,2,time_lag_for_model+1), cov_matrix(2,2,time_lag_for_model+1)];
								else
									upper_right_quadrant = cov_matrix(:,:,time_lag_for_model+1).';
								end
								
								% CCS
								phiid_atoms_CCS		= PhiIDFull_Analytical(full_time_lagged_cov, 'CCS');
								
								phiidCE_CCS(k,p)		= phiid_atoms_CCS.str + phiid_atoms_CCS.stx + ...
									phiid_atoms_CCS.sty + phiid_atoms_CCS.sts;
								phiidDC_CCS(k,p)		= phiid_atoms_CCS.str + phiid_atoms_CCS.stx + ...
									phiid_atoms_CCS.sty;
								phiidCD_CCS(k,p)		= phiidCE_CCS(k,p) - phiidDC_CCS(k,p);
								phiidDoubleRed_CCS(k,p) = phiid_atoms_CCS.rtr;
								phiidDoubleSyn_CCS(k,p) = phiid_atoms_CCS.sts;
								
								phiidUC_CCS(k,p)		= phiid_atoms_CCS.rts + phiid_atoms_CCS.xts + phiid_atoms_CCS.yts;
								phiidSyn_CCS(k,p)		= phiidUC_CCS(k,p) + phiidDC_CCS(k,p) + phiidDoubleSyn_CCS(k,p);
								phiidTransfer_CCS(k,p)	= phiid_atoms_CCS.xty + phiid_atoms_CCS.ytx;
								
								
							elseif strcmp(measures{l}, 'integrated_info_measures')
								
								%------------------------------------------------------------------------
								% INTEGRATED INFORMATION MEASURES
								%------------------------------------------------------------------------
								
								% makeLaggedCovariance calculates the covariance, but with ones in the diagonals,
								% and values ranging from -1 to 1 in the off-diagonals; this full time-lagged
								% covariance matrix will be different from the one above, as the one
								% above normalizes coupling matrices with a decay factor using specnorm(),
								% therefore, we'll use the one above to keep things consistent
								
								% full_time_lagged_cov2 = makeLaggedCovariance(coupling_matrix, noise_corr(1,2));
								
								% construct full time-lagged covariance matrix
								
								% switch off-diagonal values in upper right quadrant
								if size(cov_matrix,1) < 3
									upper_right_quadrant = [cov_matrix(1,1,time_lag_for_model+1), cov_matrix(2,1,time_lag_for_model+1); ...
										cov_matrix(1,2,time_lag_for_model+1), cov_matrix(2,2,time_lag_for_model+1)];
								else
									upper_right_quadrant = cov_matrix(:,:,time_lag_for_model+1).';
								end
								
								% put all quadrants together
								full_time_lagged_cov = [cov_matrix(:,:,1), cov_matrix(:,:,time_lag_for_model+1); upper_right_quadrant, cov_matrix(:,:,1)];
								
								for z = 1:length(calcNames);
									
									% name template to instantiate JIDT calculators
									class_template = 'infodynamics.measures.continuous.gaussian.%sCalculatorGaussian';
									
									% the following uses only the full time-lagged covariance matrix,to compute the measure
									if strcmp(calcNames{z}, 'AverageCorrelation')
										try
											diag_full_time_lagged_cov = diag(diag(full_time_lagged_cov));
											slS = (diag_full_time_lagged_cov^-0.5) * full_time_lagged_cov * (diag_full_time_lagged_cov^-0.5);
											measure_value = (sum(abs(slS(:))) - sum(diag(abs(slS))))/12.0;
											average_corr(k,p) = measure_value;
										catch
											average_corr(k,p) = NaN;
										end
									else
										try
											calc = javaObject(sprintf(class_template, calcNames{z}));
											calc.initialise(length(coupling_matrix));
											calc.setLaggedCovariance(full_time_lagged_cov);
											measure_value = calc.computeAverageLocalOfObservations();
										catch
											measure_value = NaN;
										end
										
										if strcmp(calcNames{z}, 'IntegratedInformation')
											integrated_info(k,p) = measure_value;
										elseif strcmp(calcNames{z}, 'IntegratedInteraction')
											integrated_interaction(k,p) = measure_value;
										elseif strcmp(calcNames{z}, 'DecoderIntegration')
											decoder_integration(k,p) = measure_value;
										elseif strcmp(calcNames{z}, 'CausalDensity')
											causal_density(k,p) = measure_value;
										elseif strcmp(calcNames{z}, 'IntegratedSynergy')
											integrated_synergy(k,p) = measure_value;
										elseif strcmp(calcNames{z}, 'TimeDelayedMutualInfo')
											time_delayed_mi(k,p) = measure_value;
										end
									end
								end
							end
							
						catch % ME.message
							err_count = err_count + 1;
						end
						count = count + 1;
						count_disp = ['count = ',num2str(count)];
					end
				end
			end
			
			%k_disp = ['coupling sample index = ',num2str(k)];
			%disp(k_disp);	
		end
		
		new_struct.two_norm				= new_norms;
		new_struct.spectral_radius			= spectral_radiuses;
		new_struct.coupling_matrix			= coupling_matrices;
		new_struct.noise_corr				= noise_corrs;
		new_struct.n_nodes				= n_nodes;
		new_struct.m_dim					= m_dim;
		
		if exist ('multi_info')
			new_struct.multi_info			= multi_info;
		end
		
		if exist ('phiidDoubleRed_MMI')
			new_struct.phiidDoubleRed_MMI		= phiidDoubleRed_MMI;
			new_struct.phiidDoubleSyn_MMI		= phiidDoubleSyn_MMI;
			new_struct.phiidCE_MMI			= phiidCE_MMI;
			new_struct.phiidDC_MMI			= phiidDC_MMI;
			new_struct.phiidCD_MMI			= phiidCD_MMI;
			new_struct.phiidUC_MMI			= phiidUC_MMI;
			new_struct.phiidSyn_MMI			= phiidSyn_MMI;
			new_struct.phiidTransfer_MMI		= phiidTransfer_MMI;
		end
		
		if exist ('phiidDoubleRed_CCS')
			new_struct.phiidDoubleRed_CCS		= phiidDoubleRed_CCS;
			new_struct.phiidDoubleSyn_CCS		= phiidDoubleSyn_CCS;
			new_struct.phiidCE_CCS			= phiidCE_CCS;
			new_struct.phiidDC_CCS			= phiidDC_CCS;
			new_struct.phiidCD_CCS			= phiidCD_CCS;
			new_struct.phiidUC_CCS			= phiidUC_CCS;
			new_struct.phiidSyn_CCS			= phiidSyn_CCS;
			new_struct.phiidTransfer_CCS		= phiidTransfer_CCS;
		end
		
		if exist ('average_corr')
			new_struct.average_corr			= average_corr;
		end
		
		if exist ('integrated_info')
			new_struct.integrated_info		= integrated_info;
		end
		
		if exist ('integrated_interaction')
			new_struct.integrated_interaction	= integrated_interaction;
		end
			
		if exist ('decoder_integration')
			new_struct.decoder_integration	= decoder_integration;
		end
		
		if exist ('causal_density')
			new_struct.causal_density		= causal_density;
		end
		
		if exist ('integrated_synergy')
			new_struct.integrated_synergy		= integrated_synergy;
		end
		
		if exist ('time_delayed_mi')
			new_struct.time_delayed_mi		= time_delayed_mi;
		end
		
		
		if exist('all_DD_pca') && exist('all_ShannonCE_pca') && exist('all_co_info_pca')
			new_struct.DD_pca				= all_DD_pca;
			new_struct.ShannonCE_pca		= all_ShannonCE_pca;
			new_struct.co_info_pca			= all_co_info_pca;
		end

		
		if exist('all_DD_grassmanian') && exist('all_ShannonCE_grassmanian') && exist('all_co_info_grassmanian')
			new_struct.DD_grassmanian		= all_DD_grassmanian;
			new_struct.ShannonCE_grassmanian	= all_ShannonCE_grassmanian;
			new_struct.co_info_grassmanian	= all_co_info_grassmanian;
		end
			
		results{i,j}					= new_struct;
		
		clear new_struct;
		clear new_norms;
		clear spectral_radiuses;
		clear couplings_matrices;
		clear noise_corrs;
		clear multi_info;
		
		clear phiidDoubleRed_MMI;
		clear phiidDoubleSyn_MMI;
		clear phiidCE_MMI;
		clear phiidDC_MMI;
		clear phiidCD_MMI;
		clear phiidUC_MMI;
		clear phiidSyn_MMI;
		clear phiidTransfer_MMI;
		
		clear phiidDoubleRed_CCS;
		clear phiidDoubleSyn_CCS;
		clear phiidCE_CCS;
		clear phiidDC_CCS;
		clear phiidCD_CCS;
		clear phiidUC_CCS;
		clear phiidSyn_CCS;
		clear phiidTransfer_CCS;
		
		clear all_DD_pca;
		clear all_ShannonCE_pca;
		clear all_DD_grassmanian;
		clear all_ShannonCE_grassmanian;
		clear all_co_info_pca;
		clear all_co_info_grassmanian;
		
		clear average_corr;
		clear integrated_info;
		clear integrated_interaction;
		clear decoder_integration;
		clear causal_density;
		clear integrated_synergy;
		clear time_delayed_mi;
		
		j_disp = ['RMI index = ',num2str(j)];
		disp(j_disp);
	end
		
	i_disp = ['NORM index = ',num2str(i)];
	disp(i_disp);
end
%% save into table

all_norms_str = {};
for t = 1:length(all_norms)
	all_norms_str{t} = num2str(all_norms(t));
end

all_rmi_str = {};
for e = 1:length(all_rmi)
	all_rmi_str{e} = num2str(all_rmi(e));
end

results_table = array2table(results, 'RowNames', all_norms_str, ...
	'VariableNames', all_rmi_str);

save([char(pathout_data_measures), network '_' char(filename_table) '.mat'],'results_table');
save([char(pathout_data_measures), network '_' char(filename) '.mat'],'results');

clear results