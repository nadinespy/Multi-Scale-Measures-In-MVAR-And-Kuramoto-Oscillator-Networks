%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?
% - add integrated information measures?

%% 256-NODE KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements practical synergy capacity for 256-node Kuramoto oscillators with different couplings and phase lags, two different macro variables 
% (variance of synchronies & global average pairwise synchrony between communities), and three different micro variables (thetas, cos(thetas), synchronies, and 
% binarized synchronies).

clear all;
clear java;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

pathout_data = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/'];
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/sim_time_series/'];
pathout_data_synchronies = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/synchronies/'];
pathout_data_binarized_synchronies = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/binarized_synchronies/'];
pathout_data_practical_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/practical_ce/'];
pathout_data_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/phiid_ce/'];
pathout_data_average_corr = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/average_corr/'];
pathout_data_average_cov = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/average_cov/'];

pathout_plots = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/'];
pathout_plots_practical_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/practical_ce/'];
pathout_plots_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/phiid_ce/'];
pathout_plots_sigma_met = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/sigma_met/'];
pathout_plots_sigma_chi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/sigma_chi/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/distributions/'];


%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
all_npoints = [2000]; %, 10000];
taus = [10]; %, 1];

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_random(), chimera_metastable_model())
sim_method = @chimera_metastable_model;

% network (options: '2node_mvar' for 2-node network with 100 different coupling strengths & noise correlations (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations (if choosing sim_index = 3);
% '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations (if choosing sim_index = 3));
% '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations (if choosing sim_index = 1 or 2)
% '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations (if choosing sim_index = 1 or 2)
% '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 4) 
% '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 5) 

network = '256node_kuramoto';
binarization_threshold = 0.8;

%% load files (if already existent, to, e. g., only create plots)

%{

load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');
load([pathout_data network '_all_average_corr_theta' sim_index '.mat'], 'all_average_corr_theta');
load([pathout_data network '_all_average_cov_theta' sim_index '.mat'], 'all_average_cov_theta');

synergy_capacity_practical_sigma_chi = emergence_practical.synergy_capacity_practical_sigma_chi; 
downward_causation_practical_sigma_chi = emergence_practical.downward_causation_practical_sigma_chi; 
causal_decoupling_practical_sigma_chi = emergence_practical.causal_decoupling_practical_sigma_chi;

synergy_capacity_practical_pairwise_synchrony = emergence_practical.synergy_capacity_practical_pairwise_synchrony; 
downward_causation_practical_pairwise_synchrony = emergence_practical.downward_causation_practical_pairwise_synchrony; 
causal_decoupling_practical_pairwise_synchrony = emergence_practical.causal_decoupling_practical_pairwise_synchrony;

%}

%% create coupling matrices & noise correlation vectors

% {

intra_comm_size = 32;							  % intra-community size
n_communities = 8;							    % number of communities		
	
coupling_vec = linspace(0.05, 0.9, 10);
beta_vec = linspace(0.0, 0.4, 10);				    % use beta values only up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
													% values of beta; in these cases, sigma chi will be a non-varying zero macro variable,
													% yielding erroneous values for emergence 		
	
d0 = intra_comm_size; 
d1 = intra_comm_size;							  % numbers of connections at different community levels
		
N = intra_comm_size*n_communities;			   % total number of oscillators: 8
M = n_communities;							   % number of lowest level communities (what's that?): 2
	
coupling_matrices = zeros(N,N,length(coupling_vec));
for o = 1:length(coupling_vec);
	A = coupling_vec(o);						     % was 0.2 (the higher A, the stronger the intra-community coupling strength)
	k1 = (1-A)/2;				
	k0 = 1-k1;						
	
	% build coupling matrix
	for i=1:N
		x1 = mod(ceil(i/intra_comm_size)-1,n_communities)+1;				% community number
		for j=i:N
			if i~=j	% ignore diagonals
				y1 = mod(ceil(j/intra_comm_size)-1,n_communities)+1;		% community number
				if x1 == y1		% same community
					p = d0/intra_comm_size;
					k = k0;
				else		% different communities
					p = d1/(intra_comm_size*n_communities);
					k = k1;
				end
				if rand < p
					coupling_matrices(i, j, o) = k;
					coupling_matrices(j, i, o) = k;
				end
			end
		end
	end
end

%% calculating practical measures for different macro and micro variables

% instantiate variables to store practical measures for synergistic capacity for different coupling matrices and noise correlations
nbins = 100;

rng(1);
for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for i = 1:size(coupling_matrices, 3);
		disp(i)
		coupling_string = num2str(coupling_vec(i));
		coupling_string = coupling_string(3:end);
			
		coupling_matrix = coupling_matrices(:,:,i);
			
		for j = 1:length(beta_vec)
			disp(j)

			beta = beta_vec(j);
			beta_string = num2str(beta);
			beta_string = beta_string(3:end);
				
			[thetas, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);
			save([pathout_data_sim_time_series network '_thetas_' num2str(npoints) '_' coupling_string '_' betastring'.mat'], 'thetas');
			save([pathout_data_sim_time_series network '_sigma_chi_' num2str(npoints) '_' coupling_string '_' betastring '.mat'], 'sigma_chi');
			save([pathout_data_sim_time_series network '_synchrony_' num2str(npoints) '_' coupling_string '_' betastring '.mat'], 'synchrony');
				
			for z = 1:length(taus);
				tau = taus(z);
				
				% check distributions of a subset of parameters
				if (((i == 1) && (j == 1)) || ((i == 3) && (j == 3)) || ((i == 7) && (j == 7)) || ((i == 10) && (j == 10)));
					
					r = randi([0 256],1,3);
					for k = 1:length(r);
						a_string = num2str(coupling_vec(i));
						a_string = a_string(3:end);
						histogram(thetas(k,:)', nbins);
						title(['Distribution of theta ', num2str(r(k)), ', A = ', num2str(coupling_vec(i)), ', beta = ', num2str(beta_vec(j))]);
						ylabel('frequency');
						xlabel(['values of theta ' num2str(r(k))]);
						exportgraphics(gcf, [pathout_plots_distributions network '_distr_theta_' num2str(r(k)) '_' a_string '_' num2str(npoints) '.png']);
					end
					close all;
					
					% 			figure;
					% 			plot(thetas.')
					% 			title(['Time-series of thetas, A = ', num2str(coupling_vec(i)) ', beta = ' num2str(beta_vec(j))]);
					% 			ylabel('thetas');
					% 			xlabel('Time');
					% 			exportgraphics(gcf, [pathout_plots network '_time_series_thetas_' a_string '_' sim_index '.png']);
					
				end
				
				synchronies(j,:,:) = synchrony;				% rows: betas, columns: communities, 3rd dimension: time-points
				
				% practical measures for causal emergence: variance of synchronies & global average pairwise synchrony between communities
				
				% generating global average pairwise synchrony between communities as a macro variable
				grand_mean_pairwise_synchrony = zeros(1,npoints);
				
				for h = 1:M;
					mean_pairwise_synchrony_temp = zeros(1,npoints);
					for g = 1:M;
						mean_pairwise_synchrony_temp = mean_pairwise_synchrony_temp + abs((synchrony(h,:)+synchrony(g,:))/2);
					end
					mean_pairwise_synchrony_temp = mean_pairwise_synchrony_temp/M;
					grand_mean_pairwise_synchrony = grand_mean_pairwise_synchrony + mean_pairwise_synchrony_temp;
				end
				grand_mean_pairwise_synchrony = grand_mean_pairwise_synchrony/M;
				
				% practical causal emergence with phases as micro variables
				synergy_capacity_practical_sigma_chi(i,j) = EmergencePsi(thetas', sigma_chi', tau);
				synergy_capacity_practical_pairwise_synchrony(i,j) = EmergencePsi(thetas', grand_mean_pairwise_synchrony', tau);
				downward_causation_practical_sigma_chi(i,j) = EmergenceDelta(thetas', sigma_chi', tau);
				downward_causation_practical_pairwise_synchrony(i,j) = EmergenceDelta(thetas', grand_mean_pairwise_synchrony', tau);
				causal_decoupling_practical_sigma_chi(i,j) = synergy_capacity_practical_sigma_chi(i,j) - downward_causation_practical_sigma_chi(i,j);
				causal_decoupling_practical_pairwise_synchrony(i,j) = synergy_capacity_practical_pairwise_synchrony(i,j) - downward_causation_practical_pairwise_synchrony(i,j);
				
				% practical causal emergence with synchronies as "micro" variables
				synergy_capacity_practical_sigma_chi_micro_sync(i,j) = EmergencePsi(synchrony', sigma_chi', tau);
				synergy_capacity_practical_pairwise_synchrony_micro_sync(i,j) = EmergencePsi(synchrony', grand_mean_pairwise_synchrony', tau);
				downward_causation_practical_sigma_chi_micro_sync(i,j) = EmergenceDelta(synchrony', sigma_chi', tau);
				downward_causation_practical_pairwise_synchrony_micro_sync(i,j) = EmergenceDelta(synchrony', grand_mean_pairwise_synchrony', tau);
				causal_decoupling_practical_sigma_chi_micro_sync(i,j) = synergy_capacity_practical_sigma_chi_micro_sync(i,j) - downward_causation_practical_sigma_chi_micro_sync(i,j);
				causal_decoupling_practical_pairwise_synchrony_micro_sync(i,j) = synergy_capacity_practical_pairwise_synchrony_micro_sync(i,j) - downward_causation_practical_pairwise_synchrony_micro_sync(i,j);
				
				% doing the same, albeit with binarized synchronies as the "micro" variables
				for k = 1:size(synchrony, 2);
					for l = 1:size(synchrony,1);
						if synchrony(l,k) >= binarization_threshold;
							binarized_synchrony(l,k) = 1;
						else binarized_synchrony(l,k) = 0;
						end
					end
				end
				
				binarized_synchronies(j,:,:) = binarized_synchrony;
				
				synergy_capacity_practical_sigma_chi_micro_sync_bin(i,j) = EmergencePsi(binarized_synchrony', sigma_chi', tau);
				synergy_capacity_practical_pairwise_synchrony_micro_sync_bin(i,j) = EmergencePsi(binarized_synchrony', grand_mean_pairwise_synchrony', tau);
				downward_causation_practical_sigma_chi_micro_sync_bin(i,j) = EmergenceDelta(binarized_synchrony', sigma_chi', tau);
				downward_causation_practical_pairwise_synchrony_micro_sync_bin(i,j) = EmergenceDelta(binarized_synchrony', grand_mean_pairwise_synchrony', tau);
				causal_decoupling_practical_sigma_chi_micro_sync_bin(i,j) = synergy_capacity_practical_sigma_chi_micro_sync_bin(i,j) - downward_causation_practical_sigma_chi_micro_sync_bin(i,j);
				causal_decoupling_practical_pairwise_synchrony_micro_sync_bin(i,j) = synergy_capacity_practical_pairwise_synchrony_micro_sync_bin(i,j) - downward_causation_practical_pairwise_synchrony_micro_sync_bin(i,j);
				
				% practical causal emergence with cos(phase) as the "true" micro variables (the "raw" signal)
				cos_theta = cos(thetas);
				
				synergy_capacity_practical_sigma_chi_micro_theta_cos(i,j) = EmergencePsi(cos_theta', sigma_chi', tau);
				synergy_capacity_practical_pairwise_synchrony_micro_theta_cos(i,j) = EmergencePsi(cos_theta', grand_mean_pairwise_synchrony', tau);
				downward_causation_practical_sigma_chi_micro_theta_cos(i,j) = EmergenceDelta(cos_theta', sigma_chi', tau);
				downward_causation_practical_pairwise_synchrony_micro_theta_cos(i,j) = EmergenceDelta(cos_theta', grand_mean_pairwise_synchrony', tau);
				causal_decoupling_practical_sigma_chi_micro_theta_cos(i,j) = synergy_capacity_practical_sigma_chi_micro_theta_cos(i,j) - downward_causation_practical_sigma_chi_micro_theta_cos(i,j);
				causal_decoupling_practical_pairwise_synchrony_micro_theta_cos(i,j) = synergy_capacity_practical_pairwise_synchrony_micro_theta_cos(i,j) - downward_causation_practical_pairwise_synchrony_micro_theta_cos(i,j);
				
				% average covariance/correlation matrix
				cov_theta = cov(thetas');
				all_average_cov_theta(i,j) = mean(nonzeros(tril(cov_theta,-1)), 'all');
				
				corr_theta = corrcov(cov_theta);
				all_average_corr_theta(i,j) = mean(nonzeros(tril(corr_theta,-1)), 'all');
				
			end
			
			a_string = num2str(coupling_vec(i));
			save([pathout_data_synchronies network '_synchronies_' a_string(3:end) '_' num2str(npoints) '.mat'], 'synchronies');
			save([pathout_data_binarized_synchronies network '_binarized_synchronies_' a_string(3:end) '_' num2str(npoints) '.mat'], 'binarized_synchronies');
			
		end
		
		save([pathout_data_average_corr network '_all_average_corr_theta_' num2str(npoints) '.mat'], 'all_average_corr_theta');
		save([pathout_data_average_cov network '_all_average_cov_theta_' num2str(npoints) '.mat'], 'all_average_cov_theta');
		
		%% storing practical measures for different macro variables in struct files
		
		% practical measures for different macro variables
		
		emergence_practical = [];
		
		emergence_practical.synergy_capacity_practical_sigma_chi = synergy_capacity_practical_sigma_chi;
		emergence_practical.causal_decoupling_practical_sigma_chi = causal_decoupling_practical_sigma_chi;
		emergence_practical.downward_causation_practical_sigma_chi = downward_causation_practical_sigma_chi;
		emergence_practical.synergy_capacity_practical_pairwise_synchrony = synergy_capacity_practical_pairwise_synchrony;
		emergence_practical.causal_decoupling_practical_pairwise_synchrony = causal_decoupling_practical_pairwise_synchrony;
		emergence_practical.downward_causation_practical_pairwise_synchrony = downward_causation_practical_pairwise_synchrony;
		
		emergence_practical.synergy_capacity_practical_sigma_chi_micro_sync = synergy_capacity_practical_sigma_chi_micro_sync;
		emergence_practical.causal_decoupling_practical_sigma_chi_micro_sync = causal_decoupling_practical_sigma_chi_micro_sync;
		emergence_practical.downward_causation_practical_sigma_chi_micro_sync = downward_causation_practical_sigma_chi_micro_sync;
		emergence_practical.synergy_capacity_practical_pairwise_synchrony_micro_sync = synergy_capacity_practical_pairwise_synchrony_micro_sync;
		emergence_practical.causal_decoupling_practical_pairwise_synchrony_micro_sync = causal_decoupling_practical_pairwise_synchrony_micro_sync;
		emergence_practical.downward_causation_practical_pairwise_synchrony_micro_sync = downward_causation_practical_pairwise_synchrony_micro_sync;
		
		emergence_practical.synergy_capacity_practical_sigma_chi_micro_sync_bin = synergy_capacity_practical_sigma_chi_micro_sync_bin;
		emergence_practical.causal_decoupling_practical_sigma_chi_micro_sync_bin = causal_decoupling_practical_sigma_chi_micro_sync_bin;
		emergence_practical.downward_causation_practical_sigma_chi_micro_sync_bin = downward_causation_practical_sigma_chi_micro_sync_bin;
		emergence_practical.synergy_capacity_practical_pairwise_synchrony_micro_sync_bin = synergy_capacity_practical_pairwise_synchrony_micro_sync_bin;
		emergence_practical.causal_decoupling_practical_pairwise_synchrony_micro_sync_bin = causal_decoupling_practical_pairwise_synchrony_micro_sync_bin;
		emergence_practical.downward_causation_practical_pairwise_synchrony_micro_sync_bin = downward_causation_practical_pairwise_synchrony_micro_sync_bin;
		
		emergence_practical.synergy_capacity_practical_sigma_chi_micro_theta_cos = synergy_capacity_practical_sigma_chi_micro_theta_cos;
		emergence_practical.causal_decoupling_practical_sigma_chi_micro_theta_cos = causal_decoupling_practical_sigma_chi_micro_theta_cos;
		emergence_practical.downward_causation_practical_sigma_chi_micro_theta_cos = downward_causation_practical_sigma_chi_micro_theta_cos;
		emergence_practical.synergy_capacity_practical_pairwise_synchrony_micro_theta_cos = synergy_capacity_practical_pairwise_synchrony_micro_theta_cos;
		emergence_practical.causal_decoupling_practical_pairwise_synchrony_micro_theta_cos = causal_decoupling_practical_pairwise_synchrony_micro_theta_cos;
		emergence_practical.downward_causation_practical_pairwise_synchrony_micro_theta_cos = downward_causation_practical_pairwise_synchrony_micro_theta_cos;
		
		save([pathout_data_practical_ce network '_emergence_practical_' num2str(npoints) '_' num2str(tau) '.mat'], 'emergence_practical');

	end 
	
	clear synchronies;
	clear binarized_synchronies;
	clear all_average_corr_theta;
	clear all_average_cov_theta;
	clear synergy_capacity_practical_sigma_chi;
	clear synergy_capacity_practical_pairwise_synchrony;
	clear downward_causation_practical_sigma_chi;
	clear downward_causation_practical_pairwise_synchrony;
	clear causal_decoupling_practical_sigma_chi;
	clear causal_decoupling_practical_pairwise_synchrony;
	clear synergy_capacity_practical_sigma_chi_micro_sync;
	clear synergy_capacity_practical_pairwise_synchrony_micro_sync;
	clear downward_causation_practical_sigma_chi_micro_sync;
	clear downward_causation_practical_pairwise_synchrony_micro_sync;
	clear causal_decoupling_practical_sigma_chi_micro_sync;
	clear causal_decoupling_practical_pairwise_synchrony_micro_sync;
	clear synergy_capacity_practical_sigma_chi_micro_sync_bin;
	clear synergy_capacity_practical_pairwise_synchrony_micro_sync_bin;
	clear downward_causation_practical_sigma_chi_micro_sync_bin;
	clear downward_causation_practical_pairwise_synchrony_micro_sync_bin;
	clear causal_decoupling_practical_sigma_chi_micro_sync_bin;
	clear causal_decoupling_practical_pairwise_synchrony_micro_sync_bin;
	clear synergy_capacity_practical_sigma_chi_micro_theta_cos;
	clear synergy_capacity_practical_pairwise_synchrony_micro_theta_cos;
	clear downward_causation_practical_sigma_chi_micro_theta_cos;
	clear downward_causation_practical_pairwise_synchrony_micro_theta_cos;
	clear causal_decoupling_practical_sigma_chi_micro_theta_cos;
	clear causal_decoupling_practical_pairwise_synchrony_micro_theta_cos;
	
end 
%}

%% plotting

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for z = 1:length(taus);
		tau = taus(z);
		
		load([pathout_data_practical_ce network '_emergence_practical_' num2str(npoints) '_' num2str(tau) '.mat'], 'emergence_practical');
		load([pathout_data_average_corr network '_all_average_corr_theta_' num2str(npoints) '.mat'], 'all_average_corr_theta');
		load([pathout_data_average_cov network '_all_average_cov_theta_' num2str(npoints) '.mat'], 'all_average_cov_theta');
		
		synergy_capacity_practical_sigma_chi = emergence_practical.synergy_capacity_practical_sigma_chi;
		downward_causation_practical_sigma_chi = emergence_practical.downward_causation_practical_sigma_chi;
		causal_decoupling_practical_sigma_chi = emergence_practical.causal_decoupling_practical_sigma_chi;
		synergy_capacity_practical_pairwise_synchrony = emergence_practical.synergy_capacity_practical_pairwise_synchrony;
		downward_causation_practical_pairwise_synchrony = emergence_practical.downward_causation_practical_pairwise_synchrony;
		causal_decoupling_practical_pairwise_synchrony = emergence_practical.causal_decoupling_practical_pairwise_synchrony;
		
		synergy_capacity_practical_sigma_chi_micro_sync = emergence_practical.synergy_capacity_practical_sigma_chi_micro_sync;
		causal_decoupling_practical_sigma_chi_micro_sync = emergence_practical.causal_decoupling_practical_sigma_chi_micro_sync;
		downward_causation_practical_sigma_chi_micro_sync = emergence_practical.downward_causation_practical_sigma_chi_micro_sync;
		synergy_capacity_practical_pairwise_synchrony_micro_sync = emergence_practical.synergy_capacity_practical_pairwise_synchrony_micro_sync;
		causal_decoupling_practical_pairwise_synchrony_micro_sync = emergence_practical.causal_decoupling_practical_pairwise_synchrony_micro_sync;
		downward_causation_practical_pairwise_synchrony_micro_sync = emergence_practical.downward_causation_practical_pairwise_synchrony_micro_sync;
		
		synergy_capacity_practical_sigma_chi_micro_sync_bin = emergence_practical.synergy_capacity_practical_sigma_chi_micro_sync_bin;
		causal_decoupling_practical_sigma_chi_micro_sync_bin = emergence_practical.causal_decoupling_practical_sigma_chi_micro_sync_bin;
		downward_causation_practical_sigma_chi_micro_sync_bin = emergence_practical.downward_causation_practical_sigma_chi_micro_sync_bin;
		synergy_capacity_practical_pairwise_synchrony_micro_sync_bin = emergence_practical.synergy_capacity_practical_pairwise_synchrony_micro_sync_bin;
		causal_decoupling_practical_pairwise_synchrony_micro_sync_bin = emergence_practical.causal_decoupling_practical_pairwise_synchrony_micro_sync_bin;
		downward_causation_practical_pairwise_synchrony_micro_sync_bin = emergence_practical.downward_causation_practical_pairwise_synchrony_micro_sync_bin ;
		
		synergy_capacity_practical_sigma_chi_micro_theta_cos = emergence_practical.synergy_capacity_practical_sigma_chi_micro_theta_cos;
		causal_decoupling_practical_sigma_chi_micro_theta_cos = emergence_practical.causal_decoupling_practical_sigma_chi_micro_theta_cos;
		downward_causation_practical_sigma_chi_micro_theta_cos = emergence_practical.downward_causation_practical_sigma_chi_micro_theta_cos;
		synergy_capacity_practical_pairwise_synchrony_micro_theta_cos = emergence_practical.synergy_capacity_practical_pairwise_synchrony_micro_theta_cos;
		causal_decoupling_practical_pairwise_synchrony_micro_theta_cos = emergence_practical.causal_decoupling_practical_pairwise_synchrony_micro_theta_cos;
		downward_causation_practical_pairwise_synchrony_micro_theta_cos = emergence_practical.downward_causation_practical_pairwise_synchrony_micro_theta_cos;

	
% 		% axes ticks
% 		if sim_index == '3'
% 			x_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ...
% 				'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
% 			y_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ...
% 				'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
% 		else
% 			x_axis = {'0.0', '', '0.1', '', '0.2', '', '0.3', '', '', '0.4'};
% 			y_axis = {'0.0', '', '0.24', '', '0.43', '', '0.62', '', '', '0.9'};
% 		end
		
		x_axis = {'0.0', '', '0.1', '', '0.2', '', '0.3', '', '', '0.4'};
		y_axis = {'0.0', '', '0.24', '', '0.43', '', '0.62', '', '', '0.9'};
		
		%% synergistic capacity, downward causation, causal decoupling
		
		% heatmaps using matlab built-in function:

		atoms = {synergy_capacity_practical_pairwise_synchrony, ...
			synergy_capacity_practical_sigma_chi, ...
			synergy_capacity_practical_pairwise_synchrony_micro_sync_bin, ...
			synergy_capacity_practical_sigma_chi_micro_sync_bin, ...
			synergy_capacity_practical_pairwise_synchrony_micro_sync, ...
			synergy_capacity_practical_sigma_chi_micro_sync, ...
			synergy_capacity_practical_pairwise_synchrony_micro_theta_cos, ...
			synergy_capacity_practical_sigma_chi_micro_theta_cos, ...
			downward_causation_practical_pairwise_synchrony, ...
			downward_causation_practical_sigma_chi, ...
			downward_causation_practical_pairwise_synchrony_micro_sync_bin, ...
			downward_causation_practical_sigma_chi_micro_sync_bin, ...
			downward_causation_practical_pairwise_synchrony_micro_sync, ...
			downward_causation_practical_sigma_chi_micro_sync, ...
			downward_causation_practical_pairwise_synchrony_micro_theta_cos, ...
			downward_causation_practical_sigma_chi_micro_theta_cos, ...
			causal_decoupling_practical_pairwise_synchrony, ...
			causal_decoupling_practical_sigma_chi, ...
			causal_decoupling_practical_pairwise_synchrony_micro_sync_bin, ...
			causal_decoupling_practical_sigma_chi_micro_sync_bin, ...
			causal_decoupling_practical_pairwise_synchrony_micro_sync, ...
			causal_decoupling_practical_sigma_chi_micro_sync, ...
			causal_decoupling_practical_pairwise_synchrony_micro_theta_cos, ...
			causal_decoupling_practical_sigma_chi_micro_theta_cos, ...
			all_average_cov_theta, ...
			all_average_corr_theta};
		
		file_names = {'_all_beta_coup_practical_synergy_capacity_pairwise_synchrony', ...
			'_all_beta_coup_practical_synergy_capacity_sigma_chi', ...
			'_all_beta_coup_practical_synergy_capacity_pairwise_synchrony_micro_sync_bin', ...
			'_all_beta_coup_practical_synergy_capacity_sigma_chi_micro_sync_bin', ...
			'_all_beta_coup_practical_synergy_capacity_pairwise_synchrony_micro_sync', ...
			'_all_beta_coup_practical_synergy_capacity_sigma_chi_micro_sync', ...
			'_all_beta_coup_practical_synergy_capacity_pairwise_synchrony_micro_theta_cos', ...
			'_all_beta_coup_practical_synergy_capacity_sigma_chi_micro_theta_cos', ...
			'_all_beta_coup_practical_downward_causation_pairwise_synchrony', ...
			'_all_beta_coup_practical_downward_causation_sigma_chi', ...
			'_all_beta_coup_practical_downward_causation_pairwise_synchrony_micro_sync_bin', ...
			'_all_beta_coup_practical_downward_causation_sigma_chi_micro_sync_bin', ...
			'_all_beta_coup_practical_downward_causation_pairwise_synchrony_micro_sync', ...
			'_all_beta_coup_practical_downward_causation_sigma_chi_micro_sync', ...
			'_all_beta_coup_practical_downward_causation_pairwise_synchrony_micro_theta_cos', ...
			'_all_beta_coup_practical_downward_causation_sigma_chi_micro_theta_cos', ...
			'_all_beta_coup_practical_causal_decoupling_pairwise_synchrony', ...
			'_all_beta_coup_practical_causal_decoupling_sigma_chi', ...
			'_all_beta_coup_practical_downward_causation_pairwise_synchrony_micro_sync_bin', ...
			'_all_beta_coup_practical_downward_causation_sigma_chi_micro_sync_bin', ...
			'_all_beta_coup_practical_causal_decoupling_pairwise_synchrony_micro_sync', ...
			'_all_beta_coup_practical_causal_decoupling_sigma_chi_micro_sync', ...
			'_all_beta_coup_practical_causal_decoupling_pairwise_synchrony_micro_theta_cos', ...
			'_all_beta_coup_practical_causal_decoupling_sigma_chi_micro_theta_cos', ...
			'_all_beta_coup_average_cov_theta', ...
			'_all_beta_coup_average_corr_theta'};
		
		titles = {{['synergy capacity practical pairwise synchrony micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical sigma chi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical pairwise synchrony micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical sigma chi micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical pairwise synchrony micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical sigma chi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical pairwise synchrony micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['synergy capacity practical sigma chi micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi  micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['average covariance thetas'] ['npoints = ' num2str(npoints)]}, ...
			{['average correlation thetas'] ['npoints = ' num2str(npoints)]}};


		for i = 1:size(atoms,2)
			
			figure;
			
			imagesc(atoms{i});
			colormap(bluewhitered);
			colorbar;
			
			hColorbar = colorbar;
			%set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]))
			
			xticks(1:size(x_axis, 2));
			yticks(1:size(y_axis, 2));
			
			set(gca,'TickLength',[0 0])
			yticklabels(y_axis);
			xticklabels(x_axis);
			
% 			if sim_index == '3'
% 				ylabel('zero coupling');
% 				xlabel('zero noise correlation');
% 				
% 			else
% 				ylabel('A');
% 				xlabel('beta');
% 			end
			
			ylabel('A');
			xlabel('beta');
			
			title(titles{i});
			exportgraphics(gcf, [pathout_plots network file_names{i} '_' num2str(npoints) '_' num2str(tau) '.png']);
			
		end
		
		close all;
	end 
end 

%}
		
%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 256-node kuramoto oscillators, with fixed A, and varying beta

% coupling_vec(20) = 0.2131
A = [2, 4, 6, 8, 10];

for p = 1:length(A);

g = A(p);

a_string = num2str(coupling_vec(g));
a_string = a_string(3:end);

synergy_capacity_practical_pairwise_synchrony_temp = emergence_practical.synergy_capacity_practical_pairwise_synchrony(g,:);
figure;
scatter(beta_vec, synergy_capacity_practical_pairwise_synchrony_temp, 60, 'filled');
title(['emergence capacity practical pairwise synchrony, A = ', num2str(coupling_vec(g))]);
ylabel('emergence capacity practical pairwise synchrony');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_practical_ce network '_synergy_capacity_practical_pairwise_synchrony_' a_string '.png']);

synergy_capacity_practical_sigma_chi_temp = emergence_practical.synergy_capacity_practical_sigma_chi(g,:);
figure;
scatter(beta_vec, synergy_capacity_practical_sigma_chi_temp, 60, 'filled');
title(['emergence capacity practical sigma chi, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity practical sigma chi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_practical_ce network '_synergy_capacity_practical_sigma_chi_' a_string '.png']);

load([pathout_data_practical_ce network '_synchronies_ '  a_string '.mat'], 'synchronies');
%load([PATHOUT3 network '_synchronies_ '  '21313' '_' sim_index '.mat'], 'synchronies');

% metastability (average of variance of time-series of synchronies for each beta; take whole synchrony time-series of each community,
% calculate variance for each of them (so that we'll have as many variances as communities), calculate average)
sigma_met_mean = [];
for i = 1:length(beta_vec);
	sigma_met = squeeze(synchronies(i,:,:));
	sigma_met_mean(i) = mean(var(sigma_met'));			% var() gets variance across rows
end

% chimera states (average of variance of synchronies per time-point for each beta; look at each time-points separately,
% calculate variance of synchrony across communities for each of them (so that we'll have as many variances as time-points), 
% calculate average)
sigma_chi_mean = [];
for i = 1:length(beta_vec);
	sigma_chi = squeeze(synchronies(i,:,:));
	sigma_chi_mean(i) = mean(var(sigma_chi));
end

figure;
scatter(beta_vec, sigma_chi_mean, 60, 'filled');
title(['sigma chi mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma chi mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_sigma_chi network '_sigma_chi_mean_' a_string '.png']);

figure;
scatter(beta_vec, sigma_met_mean, 60, 'filled');
title(['sigma met mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma met mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_sigma_met network '_sigma_met_mean_' a_string '.png']);

end 

close all;

%% phiid-based measures (don't work for this system size)

% some covariance matrices are not positive definite which is why some information-theoretic functions in PhiID 
% can't be computed (in this case, all terms involving t1 and t2, i.e., the indices of first & second target partition):
% h_p1t1t2 = h([p1 t1 t2]);									 
% h_p2t1t2 = h([p2 t1 t2]);
% h_t1t2 = h([t1 t2]);											
% h_p1p2t1t2 = h([p1 p2 t1 t2]);
% where h = @(idx) -log(mvnpdf(sX(idx,:)', mu(idx), S(idx,idx))); 
% h() gives the multivariate entropy, takes as an input the indices of the variables to consider in sX in PhiIDFull()

% one way to check whether a given matrix is positive definite is to see whether p is positive or not - it will be
% zero, if the matrix is positive definite, and positive, if it's not:
% [~,p] = chol(some_matrix)

% calculating information atoms

rng(1);
for i = 1:size(coupling_matrices, 3);
	
	coupling_matrix = coupling_matrices(:,:,i);
	disp(i)
	
	for j = 1:length(beta_vec)
	
		beta = beta_vec(j);
		
		[thetas, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);	
		
		% partition indices (arbitrary assigment)
		part1 = linspace(1, 128, 128);				    % indices of partition 1, e.g., in an 8-element system, [1,2,7,8]
		part2 = linspace(129, 256, 128);				  % indices of partition 2, e.g., in an 8-element system, [3,4,5,6]

		% stack data and call full PhiID function
		X1 = thetas(part1,1:end-tau);
		X2 = thetas(part2,1:end-tau);					    % Here, we define X1, X2, Y1, and Y2, so the partitions at time t (X1, X2), and the partitions at time t+1 (Y1, Y2)
            Y1 = thetas(part1,1+tau:end);
		Y2 = thetas(part2,1+tau:end);
		
		X = [X1; X2; Y1; Y2];								% stack all variables from partitions row-wise
		sX = X./repmat(std(X')', [1, npoints - tau]);
		X1 = sX(1:128,:);
		X2 = sX(129:256,:);
		Y1 = sX(257:384,:);
		Y1 = sX(385:512,:);

		% PhiID
		try
			phiid_all_beta_coup_mmi(:,i,j) = struct2array(PhiIDFull(X1, X2, Y1, Y2, 'MMI'))';
		catch 
			phiid_all_beta_coup_mmi(:,i,j) = NaN;
		end
			
		try
			phiid_all_beta_coup_ccs(:,i,j) = struct2array(PhiIDFull(X1, X2, Y1, Y2, 'ccs'))';
		catch 
			phiid_all_beta_coup_ccs(:,i,j) = NaN;
		end
			
	end
	
end 

%% storing information atoms in struct files

% information atoms: extract single atoms such that rows are the couplings, and columns the errors

% rtr:  {1}{2}-->{1}{2}
% rtx: {1}{2}-->{1}
% rty: {1}{2}-->{2}
% rts: {1}{2}-->{12}
% xtr: {1}-->{1}{2}
% xtx: {1}-->{1}
% xty: {1}-->{2} 
% xts: {1}-->{12}
% ytr: {2}-->{1}{2}
% ytx: {2}-->{1}
% yty: {2}-->{2}
% yts: {2}-->{12} 
% str: {12}-->{1}{2}
% stx: {12}-->{1} 
% sty: {12}-->{2} 
% sts: {12}-->{12}

all_atoms_beta_coup_mmi = [];
all_atoms_beta_coup_mmi.rtr = squeeze(phiid_all_beta_coup_mmi(1,:,:));
all_atoms_beta_coup_mmi.rtx = squeeze(phiid_all_beta_coup_mmi(2,:,:));
all_atoms_beta_coup_mmi.rty = squeeze(phiid_all_beta_coup_mmi(3,:,:));
all_atoms_beta_coup_mmi.rts = squeeze(phiid_all_beta_coup_mmi(4,:,:));
all_atoms_beta_coup_mmi.xtr = squeeze(phiid_all_beta_coup_mmi(5,:,:));
all_atoms_beta_coup_mmi.xtx = squeeze(phiid_all_beta_coup_mmi(6,:,:));
all_atoms_beta_coup_mmi.xty = squeeze(phiid_all_beta_coup_mmi(7,:,:));
all_atoms_beta_coup_mmi.xts = squeeze(phiid_all_beta_coup_mmi(8,:,:));
all_atoms_beta_coup_mmi.ytr = squeeze(phiid_all_beta_coup_mmi(9,:,:));
all_atoms_beta_coup_mmi.ytx = squeeze(phiid_all_beta_coup_mmi(10,:,:));
all_atoms_beta_coup_mmi.yty = squeeze(phiid_all_beta_coup_mmi(11,:,:));
all_atoms_beta_coup_mmi.yts = squeeze(phiid_all_beta_coup_mmi(12,:,:));
all_atoms_beta_coup_mmi.str = squeeze(phiid_all_beta_coup_mmi(13,:,:));
all_atoms_beta_coup_mmi.stx = squeeze(phiid_all_beta_coup_mmi(14,:,:));
all_atoms_beta_coup_mmi.sty = squeeze(phiid_all_beta_coup_mmi(15,:,:));
all_atoms_beta_coup_mmi.sts = squeeze(phiid_all_beta_coup_mmi(16,:,:));

all_atoms_beta_coup_ccs = [];
all_atoms_beta_coup_ccs.rtr = squeeze(phiid_all_beta_coup_ccs(1,:,:));
all_atoms_beta_coup_ccs.rtx = squeeze(phiid_all_beta_coup_ccs(2,:,:));
all_atoms_beta_coup_ccs.rty = squeeze(phiid_all_beta_coup_ccs(3,:,:));
all_atoms_beta_coup_ccs.rts = squeeze(phiid_all_beta_coup_ccs(4,:,:));
all_atoms_beta_coup_ccs.xtr = squeeze(phiid_all_beta_coup_ccs(5,:,:));
all_atoms_beta_coup_ccs.xtx = squeeze(phiid_all_beta_coup_ccs(6,:,:));
all_atoms_beta_coup_ccs.xty = squeeze(phiid_all_beta_coup_ccs(7,:,:));
all_atoms_beta_coup_ccs.xts = squeeze(phiid_all_beta_coup_ccs(8,:,:));
all_atoms_beta_coup_ccs.ytr = squeeze(phiid_all_beta_coup_ccs(9,:,:));
all_atoms_beta_coup_ccs.ytx = squeeze(phiid_all_beta_coup_ccs(10,:,:));
all_atoms_beta_coup_ccs.yty = squeeze(phiid_all_beta_coup_ccs(11,:,:));
all_atoms_beta_coup_ccs.yts = squeeze(phiid_all_beta_coup_ccs(12,:,:));
all_atoms_beta_coup_ccs.str = squeeze(phiid_all_beta_coup_ccs(13,:,:));
all_atoms_beta_coup_ccs.stx = squeeze(phiid_all_beta_coup_ccs(14,:,:));
all_atoms_beta_coup_ccs.sty = squeeze(phiid_all_beta_coup_ccs(15,:,:));
all_atoms_beta_coup_ccs.sts = squeeze(phiid_all_beta_coup_ccs(16,:,:));

save([pathout_data network '_all_atoms_beta_coup_ccs_' sim_index '.mat'], 'all_atoms_beta_coup_ccs');
save([pathout_data network '_all_atoms_beta_coup_mmi_' sim_index '.mat'], 'all_atoms_beta_coup_mmi');

%}

%% phiid-based synergistic/emergent capacity, downward causation, causal decoupling

% {
% calculate:
% Syn(X_t;X_t-1) (synergistic capacity of the system) 
% Un (Vt;Xt'|Xt) (causal decoupling - the top term in the lattice) 
% Un(Vt;Xt'α|Xt) (downward causation) 

% synergy (only considering the synergy that the sources have, not the target): 
% {12} --> {1}{2} +
% {12} --> {1} + 
% {12} --> {2} +
% {12} --> {12} 
 
% causal decoupling: {12} --> {12}

% downward causation: 
% {12} --> {1}{2} + 
% {12} --> {1} + 
% {12} --> {2}

% synergistic capacity
synergy_capacity_ccs = all_atoms_beta_coup_ccs.str + ...
	all_atoms_beta_coup_ccs.stx + all_atoms_beta_coup_ccs.sty + all_atoms_beta_coup_ccs.sts;

downward_causation_ccs = all_atoms_beta_coup_ccs.str + all_atoms_beta_coup_ccs.stx + all_atoms_beta_coup_ccs.sty;

synergy_capacity_mmi = all_atoms_beta_coup_mmi.str + ...
	all_atoms_beta_coup_mmi.stx + all_atoms_beta_coup_mmi.sty + all_atoms_beta_coup_mmi.sts;

downward_causation_mmi = all_atoms_beta_coup_mmi.str + all_atoms_beta_coup_mmi.stx + all_atoms_beta_coup_mmi.sty;

causal_decoupling_ccs = synergy_capacity_ccs - downward_causation_ccs;
causal_decoupling_mmi = synergy_capacity_mmi - downward_causation_mmi;

% save variables in a struct
emergence_ccs = [];
emergence_mmi = [];

emergence_ccs.synergy_capacity_ccs = synergy_capacity_ccs;
emergence_ccs.causal_decoupling_ccs = causal_decoupling_ccs;
emergence_ccs.downward_causation_ccs = downward_causation_ccs;

emergence_mmi.synergy_capacity_mmi = synergy_capacity_mmi;
emergence_mmi.causal_decoupling_mmi = causal_decoupling_mmi;
emergence_mmi.downward_causation_mmi = downward_causation_mmi;

save([pathout_data network '_emergence_ccs_' sim_index '.mat'], 'emergence_ccs');
save([pathout_data network '_emergence_mmi_' sim_index '.mat'], 'emergence_mmi');
%}

%% plotting

% axes ticks
if sim_index == '3'
	x_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ... 
		'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
	y_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ... 
		'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
else
	x_axis = {'0.0', '', '0.1', '', '0.2', '', '0.3', '', '', '0.4'};
	y_axis = {'0.0', '', '0.24', '', '0.43', '', '0.62', '', '', '0.9'};
end

% synergistic capacity, downward causation, causal decoupling

% {
% heatmaps using matlab built-in function:

atoms = {synergy_capacity_ccs, ...
	downward_causation_ccs, ...
	causal_decoupling_ccs, ...
	synergy_capacity_mmi, ...
	downward_causation_mmi, ...
	causal_decoupling_mmi};

file_names = {'_all_beta_coup_synergy_capacity_ccs_', ...
	'_all_beta_coup_downward_causation_ccs_', ...
	'_all_beta_coup_practical_causal_decoupling_ccs_', ...
	'_all_beta_coup_synergy_capacity_mmi_', ...
	'_all_beta_coup_downward_causation_mmi_', ...
	'_all_beta_coup_practical_causal_decoupling_mmi_'};

titles = {'synergy capacity ccs', ...
	'downward causation ccs', ...
	'causal decoupling ccs', ...
	'synergy capacity mmi', ...
	'downward causation mmi', ...
	'causal decoupling mmi'};

for i = 1:size(atoms,2)
	
	figure;

	imagesc(atoms{i});
	colormap(bluewhitered); 
	colorbar;
	
	hColorbar = colorbar;
	set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]))
	
	xticks(1:size(x_axis, 2));
	yticks(1:size(y_axis, 2));
	
	set(gca,'TickLength',[0 0])
	yticklabels(y_axis);
	xticklabels(x_axis);
	
	if sim_index == '3'
		ylabel('zero coupling');
		xlabel('zero noise correlation');
	
	else 
		ylabel('A');
		xlabel('beta');
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} '_' num2str(npoints) '_' num2str(tau) '.png']);

end
%}

%close all;
