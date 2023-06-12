clear all;
close all;
clear java;
clc;

cd '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/analytical/scripts'
directories = @get_mvar_directories;

n_nodes			= 2;		% number of variables in network
mdim				= 1;		% dimension of macro variable
time_lag_for_model	= 1;		% number of time-lags
n_samples_couplings	= 10;		% number of random samples for connectivity matrices
n_samples_noise_corrs	= 2;		% number of random samples for noise correlation matrices
seed				= 1;

network			= [num2str(n_nodes) 'mvar' ...	   % model name
				  '_lag' num2str(time_lag_for_model)];
			  
rng_seed(seed) % seed random number generator

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
						
						% process covariance matrix (solve DLYAP): cov_matrix(:,:,0) is the time-lagged covariance 
						% matrix between t and t-1 etc, cov_matrix(:,:,1) the time-lagged covariance 
						% matrix between t and t-2 etc.
						cov_matrix		= var_to_autocov(full_coupling_matrix,noise_corr,time_lag_for_model);	 
						total_cov(k,p)	= sum(log(diag(cov_matrix(:,:,1))));	% process total variance
						general_cov(k,p)	= log(det(cov_matrix(:,:,1)));		% process generalised variance (includes covariances)
						multi_info(k,p)	= total_cov(k,p) - general_cov(k,p);	% process multi-information (nats)				
						
						%------------------------------------------------------------------------
						% DYNAMICAL DEPENDENCE & SHANNON-BASED CAUSAL EMERGENCE
						%------------------------------------------------------------------------
						
						% run functions below without plotting and printing outputs (evalc() is a nasty solution)
						evalc('get_ssi_from_var(n_nodes, time_lag_for_model, coupling_matrix, full_coupling_matrix, spectral_radius, all_rmi(j), noise_corr)');
						evalc('get_preoptimised_dd(mdim)');
						evalc('get_optimised_dd(mdim)');
						evalc('[DD, CE] = get_ce_and_dd(mdim)');
						
						all_DD(k,p)		 = min(DD);
						all_ShannonCE(k,p) = max(CE);
						
						%------------------------------------------------------------------------
						% PHIID-BASED CAUSAL EMERGENCE, DOWNWARD CAUSATION, CAUSAL DECOUPLING
						%------------------------------------------------------------------------
						
						% construct full time-lagged covariance matrix
						
						% switch off-diagonal values in upper right quadrant 
						upper_right_quadrant = [cov_matrix(1,1,time_lag_for_model+1), cov_matrix(2,1,time_lag_for_model+1); ...
							cov_matrix(1,2,time_lag_for_model+1), cov_matrix(2,2,time_lag_for_model+1)];
						
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
						
						%------------------------------------------------------------------------
						% INTEGRATED INFORMATION MEASURES
						%------------------------------------------------------------------------
						
						calcNames = {'IntegratedInformation', 'IntegratedInteraction', 'AverageCorrelation', ...
								 'DecoderIntegration', 'CausalDensity', 'IntegratedSynergy', ...
								 'TimeDelayedMutualInfo'};
						
						% get full time-lagged correlation matrix	(makeLaggedCovariance calculates the covariance,
						% but with ones in the diagonals, and values ranging from -1 to 1 in the off-diagonals;
						% this full time-lagged covariance matrix will be different from the one above, as the one
						% above normalizes coupling matrices with a decay factor using specnorm()
						% full_time_lagged_cov2 = makeLaggedCovariance(coupling_matrix, noise_corr(1,2));
						
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
	 
					catch % ME.message
						err_count = err_count + 1;
					end
					count = count + 1;
					
					count_disp = ['count = ',num2str(count)];
				end
				
			end
		end
		
		new_struct.two_norm			= new_norms;
		new_struct.spectral_radius		= spectral_radiuses;
		new_struct.coupling_matrix		= coupling_matrices;
		new_struct.noise_corr			= noise_corrs;
		new_struct.multi_info			= multi_info;

		new_struct.phiidDoubleRed_MMI		= phiidDoubleRed_MMI;
		new_struct.phiidDoubleSyn_MMI		= phiidDoubleSyn_MMI;
		new_struct.phiidCE_MMI			= phiidCE_MMI;
		new_struct.phiidDC_MMI			= phiidDC_MMI;
		new_struct.phiidCD_MMI			= phiidCD_MMI;
		new_struct.phiidUC_MMI			= phiidUC_MMI;
		new_struct.phiidSyn_MMI			= phiidSyn_MMI;
		new_struct.phiidTransfer_MMI		= phiidTransfer_MMI;
	
		new_struct.phiidDoubleRed_CCS		= phiidDoubleRed_CCS;
		new_struct.phiidDoubleSyn_CCS		= phiidDoubleSyn_CCS;
		new_struct.phiidCE_CCS			= phiidCE_CCS;
		new_struct.phiidDC_CCS			= phiidDC_CCS;
		new_struct.phiidCD_CCS			= phiidCD_CCS;
		new_struct.phiidUC_CCS			= phiidUC_CCS;
		new_struct.phiidSyn_CCS			= phiidSyn_CCS;
		new_struct.phiidTransfer_CCS		= phiidTransfer_CCS;
		
		new_struct.DD				= all_DD;
		new_struct.ShannonCE			= all_ShannonCE;
		
		new_struct.average_corr			= average_corr;
		new_struct.integrated_info		= integrated_info;
		new_struct.integrated_interaction	= integrated_interaction;
		new_struct.decoder_integration	= decoder_integration;
		new_struct.causal_density		= causal_density;
		new_struct.integrated_synergy		= integrated_synergy;
		new_struct.time_delayed_mi		= time_delayed_mi;
		
		measures{i,j}				= new_struct;
		
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
		
		clear all_DD;
		clear all_ShannonCE; 
		
		clear average_corr;
		clear integrated_info;
		clear integrated_interaction;
		clear decoder_integration;
		clear causal_density;
		clear integrated_synergy;
		clear time_delayed_mi;
		
		j_disp = ['j = ',num2str(j)];
		disp(j_disp);
	end
	
	i_disp = ['i = ',num2str(i)];
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

measures_table = array2table(measures, 'RowNames', all_norms_str, ...
	'VariableNames', all_rmi_str);

save([char(pathout_data_measures), network '_measures_table.mat'],'measures_table');
save([char(pathout_data_measures), network '_measures.mat'],'measures');