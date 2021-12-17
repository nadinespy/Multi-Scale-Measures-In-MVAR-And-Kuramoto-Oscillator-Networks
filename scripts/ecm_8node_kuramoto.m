%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - separate calculations for phiid-based and practical measures 
% - fill struct file in a loop?
% - add integrated information measures?

%% 8-NODE KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements causal emergence for 8-node Kuramoto oscillators with different couplings and phase lags, two different macro variables 
% (variance of synchronies & global average pairwise synchrony between communities), and three different micro variables (theta, cos(theta), synchronies, and 
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
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/']; 
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/sim_time_series/'];
pathout_data_synchronies = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/synchronies/'];
pathout_data_binarized_synchronies = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/binarized_synchronies/'];
pathout_data_pract_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/practical_ce/'];
pathout_data_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/phiid_ce/'];
pathout_data_average_corr = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/average_corr/'];
pathout_data_average_cov = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/average_cov/'];

pathout_plots = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/'];
pathout_plots_pract_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/practical_ce/'];
pathout_plots_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/phiid_ce/'];
pathout_plots_sigma_met = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/sigma_met/'];
pathout_plots_sigma_chi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/sigma_chi/'];


%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
all_npoints = [2000, 10000];
taus = [1 10 100];

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_random(), chimera_metastable_model())
sim_method = @chimera_metastable_model;

% network (options: '2node_mvar' for 2-node network with 100 different coupling strengths & noise correlations (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations (if choosing sim_index = 3);
% '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations (if choosing sim_index = 3));
% '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations (if choosing sim_index = 1 or 2)
% '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations (if choosing sim_index = 1 or 2)
% '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 4) 
% '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 5) 

network = '8node_kuramoto';
binarization_threshold = 0.8;

%% load files (if already existent, to, e. g., only create plots)

%{

load([pathout_data network '_ce_ccs' sim_index '.mat'], 'ce_ccs');
load([pathout_data network '_ce_mmi' sim_index '.mat'], 'ce_mmi');
load([pathout_data network '_ce_pract ' sim_index '.mat'], 'ce_pract ');
load([pathout_data network '_all_atoms_beta_coup_ccs' sim_index '.mat'], 'all_atoms_beta_coup_ccs');
load([pathout_data network '_all_atoms_beta_coup_mmi' sim_index '.mat'], 'all_atoms_beta_coup_mmi');
load([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
load([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

synergy_capacity_mmi = ce_mmi.synergy_capacity_mmi;
dc_mmi = ce_mmi.dc_mmi;
cd_mmi = ce_mmi.cd_mmi;

synergy_capacity_ccs = ce_ccs.synergy_capacity_ccs;
dc_ccs = ce_ccs.dc_ccs;
cd_ccs = ce_ccs.cd_ccs;

synergy_capacity_pract _sigma_chi = ce_pract .synergy_capacity_pract _sigma_chi; 
dc_pract _sigma_chi = ce_pract .dc_pract _sigma_chi; 
cd_pract _sigma_chi = ce_pract .cd_pract _sigma_chi;

synergy_capacity_pract _pairwise_synchrony = ce_pract .synergy_capacity_pract _pairwise_synchrony; 
dc_pract _pairwise_synchrony = ce_pract .dc_pract _pairwise_synchrony; 
cd_pract _pairwise_synchrony = ce_pract .cd_pract _pairwise_synchrony;

%}

%% create coupling matrices & noise correlation vectors

% {

intra_comm_size = 4;							% intra-community size
n_communities = 2;							% number of communities		
	
coupling_vec = linspace(0.05, 0.9, 10);
beta_vec = linspace(0.0, 0.8, 10);				
	
d0 = intra_comm_size; 
d1 = intra_comm_size;			     % numbers of connections at different community levels
		
N = intra_comm_size*n_communities;	% total number of oscillators: 8
M = n_communities;				       % number of lowest level communities (what's that?): 2
	
coupling_matrices = zeros(N,N,length(coupling_vec));
for o = 1:length(coupling_vec);
	A = coupling_vec(o);			% was 0.2 (the higher A, the stronger the intra-community coupling strength)
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

%% calculating information atoms & practical measures

for q = 1:length(all_npoints);
	disp(q)
	npoints = all_npoints(q);
	
	for i = 1:size(coupling_matrices, 3);
		disp(i)
		coupling_string = num2str(coupling_vec(i));
		coupling_string = coupling_string(3:end);
		coupling_matrix = coupling_matrices(:,:,i);
	
		for j = 1:length(beta_vec)
			
			beta = beta_vec(j);
			beta_string = num2str(beta);
			beta_string = beta_string(3:end);
				
			[theta, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);	 
			save([pathout_data_sim_time_series network '_theta_' num2str(npoints) '_' coupling_string '_' beta_string '.mat'], 'theta');
			save([pathout_data_sim_time_series network '_sigma_chi_' num2str(npoints) '_' coupling_string '_' beta_string '.mat'], 'sigma_chi');
			save([pathout_data_sim_time_series network '_synchrony_' num2str(npoints) '_' coupling_string '_' beta_string '.mat'], 'synchrony');
			
		end 
	end 
end 


for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for z = 1:length(taus);
		tau = taus(z);

		rng(1);
		for i = 1:size(coupling_matrices, 3);
			coupling_string = num2str(coupling_vec(i));
			coupling_string = coupling_string(3:end);
			coupling_matrix = coupling_matrices(:,:,i);
	
			for j = 1:length(beta_vec)
				disp(j)
			
				beta = beta_vec(j);
				beta_string = num2str(beta);
				beta_string = beta_string(3:end);
				
				load([pathout_data_sim_time_series network '_theta_' num2str(npoints) '_' coupling_string '_' beta_string '.mat'], 'theta');
				load([pathout_data_sim_time_series network '_sigma_chi_' num2str(npoints) '_' coupling_string '_' beta_string '.mat'], 'sigma_chi');
				load([pathout_data_sim_time_series network '_synchrony_' num2str(npoints) '_' coupling_string '_' beta_string '.mat'], 'synchrony');
			
				% PhiID
				% thetas
				try
					phiid_all_beta_coup_mmi_theta(:,i,j) = struct2array(PhiIDFull(theta, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_theta(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_theta(:,i,j) = struct2array(PhiIDFull(theta, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_theta(:,i,j) = NaN;
				end
				
				% cos(theta)
				cos_theta = cos(theta);
				try
					phiid_all_beta_coup_mmi_cos_theta(:,i,j) = struct2array(PhiIDFull(cos_theta, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_cos_theta(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_cos_theta(:,i,j) = struct2array(PhiIDFull(cos_theta, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_cos_theta(:,i,j) = NaN;
				end
				
				% synchronies
				try
					phiid_all_beta_coup_mmi_sync(:,i,j) = struct2array(PhiIDFull(synchrony, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_sync(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_sync(:,i,j) = struct2array(PhiIDFull(synchrony, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_sync(:,i,j) = NaN;
				end
				
				% binarized synchronies
				for k = 1:size(synchrony, 2);
					for l = 1:size(synchrony,1);
						if synchrony(l,k) >= binarization_threshold;
							binarized_synchrony(l,k) = 1;
						else binarized_synchrony(l,k) = 0;
						end
					end
				end
				
				bin_synchronies(j,:,:) = binarized_synchrony;
				
				try
					phiid_all_beta_coup_mmi_bin_sync(:,i,j) = struct2array(PhiIDFull(binarized_synchrony, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_bin_sync(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_bin_sync(:,i,j) = struct2array(PhiIDFull(binarized_synchrony, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_bin_sync(:,i,j) = NaN;
				end
				
				synchronies(j,:,:) = synchrony;				% rows: betas, columns: communities, 3rd dimension: time-points
				
				% practical measures for causal emergence: variance of synchronies & global average pairwise synchrony between communities
				grand_mean_pair_sync = zeros(1,npoints);
				
				for h = 1:M;
					mean_pair_sync_temp = zeros(1,npoints);
					for g = 1:M;
						mean_pair_sync_temp = mean_pair_sync_temp + abs((synchrony(h,:)+synchrony(g,:))/2);
					end
					mean_pair_sync_temp = mean_pair_sync_temp/M;
					grand_mean_pair_sync = grand_mean_pair_sync + mean_pair_sync_temp;
				end
				grand_mean_pair_sync = grand_mean_pair_sync/M;
				
				% practical causal emergence with thetas as micro variables
				ce_pract_sigma_chi_micro_theta(i,j) = EmergencePsi(theta', sigma_chi');
				ce_pract_pair_sync_micro_theta(i,j) = EmergencePsi(theta', grand_mean_pair_sync');
				dc_pract_sigma_chi_micro_theta(i,j) = EmergenceDelta(theta', sigma_chi');
				dc_pract_pair_sync_micro_theta(i,j) = EmergenceDelta(theta', grand_mean_pair_sync');
				cd_pract_sigma_chi_micro_theta(i,j) = ce_pract_sigma_chi_micro_theta(i,j) - dc_pract_sigma_chi_micro_theta(i,j);
				cd_pract_pair_sync_micro_theta(i,j) = ce_pract_pair_sync_micro_theta(i,j) - dc_pract_pair_sync_micro_theta(i,j);
				
				% practical causal emergence with synchronies as "micro" variables
				ce_pract_sigma_chi_micro_sync(i,j) = EmergencePsi(synchrony', sigma_chi', tau);
				ce_pract_pair_sync_micro_sync(i,j) = EmergencePsi(synchrony', grand_mean_pair_sync', tau);
				dc_pract_sigma_chi_micro_sync(i,j) = EmergenceDelta(synchrony', sigma_chi', tau);
				dc_pract_pair_sync_micro_sync(i,j) = EmergenceDelta(synchrony', grand_mean_pair_sync', tau);
				cd_pract_sigma_chi_micro_sync(i,j) = ce_pract_sigma_chi_micro_sync(i,j) - dc_pract_sigma_chi_micro_sync(i,j);
				cd_pract_pair_sync_micro_sync(i,j) = ce_pract_pair_sync_micro_sync(i,j) - dc_pract_pair_sync_micro_sync(i,j);
				
				% doing the same, albeit with binarized synchronies as the "micro" variables
				ce_pract_sigma_chi_micro_bin_sync(i,j) = EmergencePsi(binarized_synchrony', sigma_chi', tau);
				ce_pract_pair_sync_micro_bin_sync(i,j) = EmergencePsi(binarized_synchrony', grand_mean_pair_sync', tau);
				dc_pract_sigma_chi_micro_bin_sync(i,j) = EmergenceDelta(binarized_synchrony', sigma_chi', tau);
				dc_pract_pair_sync_micro_bin_sync(i,j) = EmergenceDelta(binarized_synchrony', grand_mean_pair_sync', tau);
				cd_pract_sigma_chi_micro_bin_sync(i,j) = ce_pract_sigma_chi_micro_bin_sync(i,j) - dc_pract_sigma_chi_micro_bin_sync(i,j);
				cd_pract_pair_sync_micro_bin_sync(i,j) = ce_pract_pair_sync_micro_bin_sync(i,j) - dc_pract_pair_sync_micro_bin_sync(i,j);
				
				% practical causal emergence with cos(phase) as the "true" micro variables (the "raw" signal)	
				ce_pract_sigma_chi_micro_cos_theta(i,j) = EmergencePsi(cos_theta', sigma_chi', tau);
				ce_pract_pair_sync_micro_cos_theta(i,j) = EmergencePsi(cos_theta', grand_mean_pair_sync', tau);
				dc_pract_sigma_chi_micro_cos_theta(i,j) = EmergenceDelta(cos_theta', sigma_chi', tau);
				dc_pract_pair_sync_micro_cos_theta(i,j) = EmergenceDelta(cos_theta', grand_mean_pair_sync', tau);
				cd_pract_sigma_chi_micro_cos_theta(i,j) = ce_pract_sigma_chi_micro_cos_theta(i,j) - dc_pract_sigma_chi_micro_cos_theta(i,j);
				cd_pract_pair_sync_micro_cos_theta(i,j) = ce_pract_pair_sync_micro_cos_theta(i,j) - dc_pract_pair_sync_micro_cos_theta(i,j);
				
				% average covariance/correlation matrix
				cov_theta = cov(theta');
				all_average_cov_theta(i,j) = mean(nonzeros(tril(cov_theta,-1)), 'all');
				
				corr_theta = corrcov(cov_theta);
				all_average_corr_theta(i,j) = mean(nonzeros(tril(corr_theta,-1)), 'all');
			
			end
			
			a_string = num2str(coupling_vec(i));
			save([pathout_data_synchronies network '_synchronies_'  a_string(3:end) '_' num2str(npoints) '.mat'], 'synchronies');
			save([pathout_data_binarized_synchronies network '_bin_synchronies_' a_string(3:end) '_' num2str(npoints) '.mat'], 'bin_synchronies');
			
		end

		save([pathout_data_average_corr network '_all_average_corr_theta_' num2str(npoints) '.mat'], 'all_average_corr_theta');
		save([pathout_data_average_cov network '_all_average_cov_theta_' num2str(npoints) '.mat'], 'all_average_cov_theta');

		%% storing information atoms & practical measures for different macro variables in struct files

		% practical measures for different macro variables

		ce_practical = [];
		
		ce_practical.ce_pract_sigma_chi_micro_theta= ce_pract_sigma_chi_micro_theta;
		ce_practical.cd_pract_sigma_chi_micro_theta = cd_pract_sigma_chi_micro_theta;
		ce_practical.dc_pract_sigma_chi_micro_theta = dc_pract_sigma_chi_micro_theta;
		ce_practical.ce_pract_pair_sync_micro_theta = ce_pract_pair_sync_micro_theta;
		ce_practical.cd_pract_pair_sync_micro_theta = cd_pract_pair_sync_micro_theta;
		ce_practical.dc_pract_pair_sync_micro_theta = dc_pract_pair_sync_micro_theta;
		
		ce_practical.ce_pract_sigma_chi_micro_sync = ce_pract_sigma_chi_micro_sync;
		ce_practical.cd_pract_sigma_chi_micro_sync = cd_pract_sigma_chi_micro_sync;
		ce_practical.dc_pract_sigma_chi_micro_sync = dc_pract_sigma_chi_micro_sync;
		ce_practical.ce_pract_pair_sync_micro_sync = ce_pract_pair_sync_micro_sync;
		ce_practical.cd_pract_pair_sync_micro_sync = cd_pract_pair_sync_micro_sync;
		ce_practical.dc_pract_pair_sync_micro_sync = dc_pract_pair_sync_micro_sync;
		
		ce_practical.ce_pract_sigma_chi_micro_bin_sync = ce_pract_sigma_chi_micro_bin_sync;
		ce_practical.cd_pract_sigma_chi_micro_bin_sync = cd_pract_sigma_chi_micro_bin_sync;
		ce_practical.dc_pract_sigma_chi_micro_bin_sync = dc_pract_sigma_chi_micro_bin_sync;
		ce_practical.ce_pract_pair_sync_micro_bin_sync = ce_pract_pair_sync_micro_bin_sync;
		ce_practical.cd_pract_pair_sync_micro_bin_sync = cd_pract_pair_sync_micro_bin_sync;
		ce_practical.dc_pract_pair_sync_micro_bin_sync = dc_pract_pair_sync_micro_bin_sync;
		
		ce_practical.ce_pract_sigma_chi_micro_cos_theta = ce_pract_sigma_chi_micro_cos_theta;
		ce_practical.cd_pract_sigma_chi_micro_cos_theta = cd_pract_sigma_chi_micro_cos_theta;
		ce_practical.dc_pract_sigma_chi_micro_cos_theta = dc_pract_sigma_chi_micro_cos_theta;
		ce_practical.ce_pract_pair_sync_micro_cos_theta = ce_pract_pair_sync_micro_cos_theta;
		ce_practical.cd_pract_pair_sync_micro_cos_theta = cd_pract_pair_sync_micro_cos_theta;
		ce_practical.dc_pract_pair_sync_micro_cos_theta = dc_pract_pair_sync_micro_cos_theta;
		
		save([pathout_data_pract_ce network '_ce_practical_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_practical');
		
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
		
		[all_atoms_beta_coup_mmi_theta, all_atoms_beta_coup_ccs_theta] =  store_atoms_in_struct(phiid_all_beta_coup_mmi_theta, phiid_all_beta_coup_ccs_theta);
		[all_atoms_beta_coup_mmi_cos_theta, all_atoms_beta_coup_ccs_cos_theta] =  store_atoms_in_struct(phiid_all_beta_coup_mmi_cos_theta, phiid_all_beta_coup_ccs_cos_theta);
		[all_atoms_beta_coup_mmi_sync, all_atoms_beta_coup_ccs_sync] =  store_atoms_in_struct(phiid_all_beta_coup_mmi_sync, phiid_all_beta_coup_ccs_sync);
		[all_atoms_beta_coup_mmi_bin_sync, all_atoms_beta_coup_ccs_bin_sync] =  store_atoms_in_struct(phiid_all_beta_coup_mmi_bin_sync, phiid_all_beta_coup_ccs_bin_sync);
		
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_theta');
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_theta');
		
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_cos_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_cos_theta');
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_cos_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_cos_theta');

		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_sync');
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_sync');
		
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_bin_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_bin_sync');
		save([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_bin_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_bin_sync');

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
		
		% causal emergence

		[ce_mmi_theta, ...
		dc_mmi_theta, ...
		cd_mmi_theta, ...
		ce_ccs_theta, ...
		dc_ccs_theta, ...
		cd_ccs_theta] = phiid_causal_emergence(all_atoms_beta_coup_mmi_theta, all_atoms_beta_coup_ccs_theta);

		[ce_mmi_cos_theta, ...
		dc_mmi_cos_theta, ...
		cd_mmi_cos_theta, ...
		ce_ccs_cos_theta, ...
		dc_ccs_cos_theta, ...
		cd_ccs_cos_theta] = phiid_causal_emergence(all_atoms_beta_coup_mmi_cos_theta, all_atoms_beta_coup_ccs_cos_theta);

		[ce_mmi_sync, ...
		dc_mmi_sync, ...
		cd_mmi_sync, ...
		ce_ccs_sync, ...
		dc_ccs_sync, ...
		cd_ccs_sync] = phiid_causal_emergence(all_atoms_beta_coup_mmi_sync, all_atoms_beta_coup_ccs_sync);

		[ce_mmi_bin_sync, ...
		dc_mmi_bin_sync, ...
		cd_mmi_bin_sync, ...
		ce_ccs_bin_sync, ...
		dc_ccs_bin_sync, ...
		cd_ccs_bin_sync] = phiid_causal_emergence(all_atoms_beta_coup_mmi_bin_sync, all_atoms_beta_coup_ccs_bin_sync);
		
		% save variables in a struct

		[ce_ccs, ce_mmi, ce_practical] = store_emergence_in_struct2(...
			ce_mmi_theta, ce_mmi_cos_theta, ce_mmi_sync, ce_mmi_bin_sync,...
			cd_mmi_theta, cd_mmi_cos_theta, cd_mmi_sync, cd_mmi_bin_sync,...
			dc_mmi_theta, dc_mmi_cos_theta, dc_mmi_sync, dc_mmi_bin_sync, ...
			ce_ccs_theta, ce_ccs_cos_theta, ce_ccs_sync, ce_ccs_bin_sync,...
			cd_ccs_theta, cd_ccs_cos_theta, cd_ccs_sync, cd_ccs_bin_sync,...
			dc_ccs_theta, dc_ccs_cos_theta, dc_ccs_sync, dc_ccs_bin_sync, ...
			ce_pract_sigma_chi_micro_theta, ce_pract_sigma_chi_micro_cos_theta, ce_pract_sigma_chi_micro_sync, ce_pract_sigma_chi_micro_bin_sync,...
			cd_pract_sigma_chi_micro_theta, cd_pract_sigma_chi_micro_cos_theta, cd_pract_sigma_chi_micro_sync, cd_pract_sigma_chi_micro_bin_sync,...
			dc_pract_sigma_chi_micro_theta, dc_pract_sigma_chi_micro_cos_theta, dc_pract_sigma_chi_micro_sync, dc_pract_sigma_chi_micro_bin_sync, ...
			ce_pract_pair_sync_micro_theta, ce_pract_pair_sync_micro_cos_theta, ce_pract_pair_sync_micro_sync, ce_pract_pair_sync_micro_bin_sync,...
			cd_pract_pair_sync_micro_theta, cd_pract_pair_sync_micro_cos_theta, cd_pract_pair_sync_micro_sync, cd_pract_pair_sync_micro_bin_sync,...
			dc_pract_pair_sync_micro_theta, dc_pract_pair_sync_micro_cos_theta, dc_pract_pair_sync_micro_sync, dc_pract_pair_sync_micro_bin_sync);
		
		save([pathout_data_phiid_ce network '_ce_ccs_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_ccs');
		save([pathout_data_phiid_ce network '_ce_mmi_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_mmi');
		save([pathout_data_pract_ce network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_practical');

	end 
	
	clear ce_ccs;
	clear ce_mmi;
	clear ce_practical;
	clear all_atoms_beta_coup_ccs_theta;
	clear all_atoms_beta_coup_mmi_theta;
	clear all_atoms_beta_coup_ccs_cos_theta;
	clear all_atoms_beta_coup_mmi_cos_theta;
	clear all_atoms_beta_coup_ccs_sync;
	clear all_atoms_beta_coup_mmi_sync;
	clear all_atoms_beta_coup_ccs_bin_sync;
	clear all_atoms_beta_coup_mmi_bin_sync;
	clear synchrony;
	clear binarized_synchrony;
	clear synchronies;
	clear bin_synchronies;
	clear all_average_corr_theta;
	clear all_average_cov_theta;
	clear ce_pract_sigma_chi_micro_theta;
	clear ce_pract_pair_sync_micro_theta;
	clear dc_pract_sigma_chi_micro_theta;
	clear dc_pract_pair_sync_micro_theta;
	clear cd_pract_sigma_chi_micro_theta;
	clear cd_pract_pair_sync_micro_theta;
	clear ce_pract_sigma_chi_micro_sync;
	clear ce_pract_pair_sync_micro_sync;
	clear dc_pract_sigma_chi_micro_sync;
	clear dc_pract_pair_sync_micro_sync;
	clear cd_pract_sigma_chi_micro_sync;
	clear cd_pract_pair_sync_micro_sync;
	clear ce_pract_sigma_chi_micro_bin_sync;
	clear ce_pract_pair_sync_micro_bin_sync;
	clear dc_pract_sigma_chi_micro_bin_sync;
	clear dc_pract_pair_sync_micro_bin_sync;
	clear cd_pract_sigma_chi_micro_bin_sync;
	clear cd_pract_pair_sync_micro_bin_sync;
	clear ce_pract_sigma_chi_micro_cos_theta;
	clear ce_pract_pair_sync_micro_cos_theta;
	clear dc_pract_sigma_chi_micro_cos_theta;
	clear dc_pract_pair_sync_micro_cos_theta;
	clear cd_pract_sigma_chi_micro_cos_theta;
	clear cd_pract_pair_sync_micro_cos_theta;
	
end 

%}


%% plotting

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for z = 1:length(taus);
		tau = taus(z);
		
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_theta');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_theta');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_cos_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_cos_theta');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_cos_theta_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_cos_theta');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_sync');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_sync');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_ccs_bin_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_ccs_bin_sync');
		load([pathout_data_phiid_ce network '_all_atoms_beta_coup_mmi_bin_sync_' num2str(npoints) '_' num2str(tau) '.mat'], 'all_atoms_beta_coup_mmi_bin_sync');
		
		load([pathout_data_phiid_ce network '_ce_ccs_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_ccs');
		load([pathout_data_phiid_ce network '_ce_mmi_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_mmi');
		load([pathout_data_pract_ce network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], 'ce_practical');
		load([pathout_data_average_corr network '_all_average_corr_theta_' num2str(npoints) '.mat'], 'all_average_corr_theta');
		load([pathout_data_average_cov network '_all_average_cov_theta_' num2str(npoints) '.mat'], 'all_average_cov_theta');
		
		
		ce_ccs.ce_ccs_theta = ce_ccs_theta;
		ce_ccs.cd_ccs_theta = cd_ccs_theta;
		ce_ccs.dc_ccs_theta = dc_ccs_theta;
		ce_mmi.ce_mmi_theta = ce_mmi_theta;
		ce_mmi.cd_mmi_theta = cd_mmi_theta;
		ce_mmi.dc_mmi_theta = dc_mmi_theta;
		
		ce_ccs.ce_ccs_cos_theta = ce_ccs_cos_theta;
		ce_ccs.cd_ccs_cos_theta = cd_ccs_cos_theta;
		ce_ccs.dc_ccs_cos_theta = dc_ccs_cos_theta;
		ce_mmi.ce_mmi_cos_theta = ce_mmi_cos_theta;
		ce_mmi.cd_mmi_cos_theta = cd_mmi_cos_theta;
		ce_mmi.dc_mmi_cos_theta = dc_mmi_cos_theta;
		
		ce_ccs.ce_ccs_sync = ce_ccs_sync;
		ce_ccs.cd_ccs_sync = cd_ccs_sync;
		ce_ccs.dc_ccs_sync = dc_ccs_sync;
		ce_mmi.ce_mmi_sync = ce_mmi_sync;
		ce_mmi.cd_mmi_sync = cd_mmi_sync;
		ce_mmi.dc_mmi_sync = dc_mmi_sync;
		
		ce_ccs.ce_ccs_bin_sync = ce_ccs_bin_sync;
		ce_ccs.cd_ccs_bin_sync = cd_ccs_bin_sync;
		ce_ccs.dc_ccs_bin_sync = dc_ccs_bin_sync;
		ce_mmi.ce_mmi_bin_sync = ce_mmi_bin_sync;
		ce_mmi.cd_mmi_bin_sync = cd_mmi_bin_sync;
		ce_mmi.dc_mmi_bin_sync = dc_mmi_bin_sync;
		

		ce_pract_sigma_chi_micro_theta= ce_practical.ce_pract_sigma_chi_micro_theta;
		dc_pract_sigma_chi_micro_theta = ce_practical.dc_pract_sigma_chi_micro_theta;
		cd_pract_sigma_chi_micro_theta = ce_practical.cd_pract_sigma_chi_micro_theta;
		ce_pract_pair_sync_micro_theta = ce_practical.ce_pract_pair_sync_micro_theta;
		dc_pract_pair_sync_micro_theta = ce_practical.dc_pract_pair_sync_micro_theta;
		cd_pract_pair_sync_micro_theta = ce_practical.cd_pract_pair_sync_micro_theta;
		
		ce_pract_sigma_chi_micro_sync = ce_practical.ce_pract_sigma_chi_micro_sync;
		cd_pract_sigma_chi_micro_sync = ce_practical.cd_pract_sigma_chi_micro_sync;
		dc_pract_sigma_chi_micro_sync = ce_practical.dc_pract_sigma_chi_micro_sync;
		ce_pract_pair_sync_micro_sync = ce_practical.ce_pract_pair_sync_micro_sync;
		cd_pract_pair_sync_micro_sync = ce_practical.cd_pract_pair_sync_micro_sync;
		dc_pract_pair_sync_micro_sync = ce_practical.dc_pract_pair_sync_micro_sync;
		
		ce_pract_sigma_chi_micro_bin_sync = ce_practical.ce_pract_sigma_chi_micro_bin_sync;
		cd_pract_sigma_chi_micro_bin_sync = ce_practical.cd_pract_sigma_chi_micro_bin_sync;
		dc_pract_sigma_chi_micro_bin_sync = ce_practical.dc_pract_sigma_chi_micro_bin_sync;
		ce_pract_pair_sync_micro_bin_sync = ce_practical.ce_pract_pair_sync_micro_bin_sync;
		cd_pract_pair_sync_micro_bin_sync = ce_practical.cd_pract_pair_sync_micro_bin_sync;
		dc_pract_pair_sync_micro_bin_sync = ce_practical.dc_pract_pair_sync_micro_bin_sync ;
		
		ce_pract_sigma_chi_micro_cos_theta = ce_practical.ce_pract_sigma_chi_micro_cos_theta;
		cd_pract_sigma_chi_micro_cos_theta = ce_practical.cd_pract_sigma_chi_micro_cos_theta;
		dc_pract_sigma_chi_micro_cos_theta = ce_practical.dc_pract_sigma_chi_micro_cos_theta;
		ce_pract_pair_sync_micro_cos_theta = ce_practical.ce_pract_pair_sync_micro_cos_theta;
		cd_pract_pair_sync_micro_cos_theta = ce_practical.cd_pract_pair_sync_micro_cos_theta;
		dc_pract_pair_sync_micro_cos_theta = ce_practical.dc_pract_pair_sync_micro_cos_theta;
		
		% axes ticks
		x_axis = {'0.0', '', '0.1', '', '0.2', '', '0.3', '', '', '0.4'};
		y_axis = {'0.0', '', '0.24', '', '0.43', '', '0.62', '', '', '0.9'};


		%% double-redundancy & double-synergy
		
		% {
		% heatmaps using imagesc()
		
		atoms = {all_atoms_beta_coup_ccs_theta.rtr, ...
			all_atoms_beta_coup_mmi_theta.rtr, ...
			all_atoms_beta_coup_ccs_cos_theta.rtr, ...
			all_atoms_beta_coup_mmi_cos_theta.rtr, ...
			all_atoms_beta_coup_ccs_sync.rtr, ...
			all_atoms_beta_coup_mmi_sync.rtr, ...
			all_atoms_beta_coup_ccs_bin_sync.rtr, ...
			all_atoms_beta_coup_mmi_bin_sync.rtr};
		
		file_names = {'_all_beta_coup_ccs_rtr_theta', ...
			'_all_beta_coup_mmi_rtr_theta', ...
			'_all_beta_coup_ccs_rtr_cos_theta', ...
			'_all_beta_coup_mmi_rtr_cos_theta', ...
			'_all_beta_coup_ccs_rtr_sync', ...
			'_all_beta_coup_mmi_rtr_sync', ...
			'_all_beta_coup_ccs_rtr_bin_sync', ...
			'_all_beta_coup_mmi_rtr_bin_sync'};
		
		titles = {{['double-redundancy ccs micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy ccs micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy ccs micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy ccs micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		% cmin = min([min(atoms{1}(:)), min(atoms{2}(:)), min(atoms{3}(:)), min(atoms{4}(:))]);
		% cmax = max([max(atoms{1}(:)), max(atoms{2}(:)), max(atoms{3}(:)), max(atoms{4}(:))]);

		for i = 1:size(atoms,2)
			
			figure;
			
			imagesc(atoms{i});
			colormap(bluewhitered);
			colorbar;
			%hColorbar = colorbar;
			%set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]))
			
			xticks(1:size(x_axis, 2));
			yticks(1:size(y_axis, 2));
			
			set(gca,'TickLength',[0 0])
			yticklabels(y_axis);
			xticklabels(x_axis);
			
			ylabel('A');
			xlabel('beta');
			
			title(titles{i});
			exportgraphics(gcf, [pathout_plots network file_names{i} '.png']);
			
		end
		
		close all

		%}
		
		%% synergistic capacity, downward causation, causal decoupling
		
		% {
		% heatmaps using matlab built-in function:
		
		atoms = {ce_ccs_theta, ...
			ce_mmi_theta, ...
			ce_ccs_cos_theta, ...
			ce_mmi_cos_theta, ...
			ce_ccs_sync, ...
			ce_mmi_sync, ...
			ce_ccs_bin_sync, ...
			ce_mmi_bin_sync, ...
			ce_pract_pair_sync_micro_theta, ...
			ce_pract_sigma_chi_micro_theta, ...
			ce_pract_pair_sync_micro_bin_sync, ...
			ce_pract_sigma_chi_micro_bin_sync, ...
			ce_pract_pair_sync_micro_sync, ...
			ce_pract_sigma_chi_micro_sync, ...
			ce_pract_pair_sync_micro_cos_theta, ...
			ce_pract_sigma_chi_micro_cos_theta, ...
			dc_ccs_theta, ...
			dc_mmi_theta, ...
			dc_ccs_cos_theta, ...
			dc_mmi_cos_theta, ...
			dc_ccs_sync, ...
			dc_mmi_sync, ...
			dc_ccs_bin_sync, ...
			dc_mmi_bin_sync, ...
			dc_pract_pair_sync_micro_theta, ...
			dc_pract_sigma_chi_micro_theta, ...
			dc_pract_pair_sync_micro_bin_sync, ...
			dc_pract_sigma_chi_micro_bin_sync, ...
			dc_pract_pair_sync_micro_sync, ...
			dc_pract_sigma_chi_micro_sync, ...
			dc_pract_pair_sync_micro_cos_theta, ...
			dc_pract_sigma_chi_micro_cos_theta, ...
			cd_ccs_theta, ...
			cd_mmi_theta, ...
			cd_ccs_cos_theta, ...
			cd_mmi_cos_theta, ...
			cd_ccs_sync, ...
			cd_mmi_sync, ...
			cd_ccs_bin_sync, ...
			cd_mmi_bin_sync, ...
			cd_pract_pair_sync_micro_theta, ...
			cd_pract_sigma_chi_micro_theta, ...
			cd_pract_pair_sync_micro_bin_sync, ...
			cd_pract_sigma_chi_micro_bin_sync, ...
			cd_pract_pair_sync_micro_sync, ...
			cd_pract_sigma_chi_micro_sync, ...
			cd_pract_pair_sync_micro_cos_theta, ...
			cd_pract_sigma_chi_micro_cos_theta, ...
			all_average_cov_theta, ...
			all_average_corr_theta};
		
		file_names = {'_all_beta_coup_ccs_ce_theta', ...
			'_all_beta_coup_mmi_ce_theta', ...
			'_all_beta_coup_ccs_ce_cos_theta', ...
			'_all_beta_coup_mmi_ce_cos_theta', ...
			'_all_beta_coup_ccs_ce_sync', ...
			'_all_beta_coup_mmi_ce_sync', ...
			'_all_beta_coup_ccs_ce_bin_sync', ...
			'_all_beta_coup_mmi_ce_bin_sync', ...
			'_all_beta_coup_pract_ce_pair_sync_micro_theta', ...
			'_all_beta_coup_pract_ce_sigma_chi_micro_theta', ...
			'_all_beta_coup_pract_ce_pair_sync_micro_bin_sync', ...
			'_all_beta_coup_pract_ce_sigma_chi_micro_bin_sync', ...
			'_all_beta_coup_pract_ce_pair_sync_micro_sync', ...
			'_all_beta_coup_pract_ce_sigma_chi_micro_sync', ...
			'_all_beta_coup_pract_ce_pair_sync_micro_cos_theta', ...
			'_all_beta_coup_pract_ce_sigma_chi_micro_cos_theta', ...
			'_all_beta_coup_ccs_dc_theta', ...
			'_all_beta_coup_mmi_dc_theta', ...
			'_all_beta_coup_ccs_dc_cos_theta', ...
			'_all_beta_coup_mmi_dc_cos_theta', ...
			'_all_beta_coup_ccs_dc_sync', ...
			'_all_beta_coup_mmi_dc_sync', ...
			'_all_beta_coup_ccs_dc_bin_sync', ...
			'_all_beta_coup_mmi_dc_bin_sync', ...
			'_all_beta_coup_pract_dc_pair_sync_micro_theta', ...
			'_all_beta_coup_pract_dc_sigma_chi_micro_theta', ...
			'_all_beta_coup_pract_dc_pair_sync_micro_bin_sync', ...
			'_all_beta_coup_pract_dc_sigma_chi_micro_bin_sync', ...
			'_all_beta_coup_pract_dc_pair_sync_micro_sync', ...
			'_all_beta_coup_pract_dc_sigma_chi_micro_sync', ...
			'_all_beta_coup_pract_dc_pair_sync_micro_cos_theta', ...
			'_all_beta_coup_pract_dc_sigma_chi_micro_cos_theta', ...
			'_all_beta_coup_ccs_cd_theta', ...
			'_all_beta_coup_mmi_cd_theta', ...
			'_all_beta_coup_ccs_cd_cos_theta', ...
			'_all_beta_coup_mmi_cd_cos_theta', ...
			'_all_beta_coup_ccs_cd_sync', ...
			'_all_beta_coup_mmi_cd_sync', ...
			'_all_beta_coup_ccs_cd_bin_sync', ...
			'_all_beta_coup_mmi_cd_bin_sync', ...
			'_all_beta_coup_pract_cd_pair_sync_micro_theta', ...
			'_all_beta_coup_pract_cd_sigma_chi_micro_theta', ...
			'_all_beta_coup_pract_dc_pair_sync_micro_bin_sync', ...
			'_all_beta_coup_pract_dc_sigma_chi_micro_bin_sync', ...
			'_all_beta_coup_pract_cd_pair_sync_micro_sync', ...
			'_all_beta_coup_pract_cd_sigma_chi_micro_sync', ...
			'_all_beta_coup_pract_cd_pair_sync_micro_cos_theta', ...
			'_all_beta_coup_pract_cd_sigma_chi_micro_cos_theta', ...
			'_all_beta_coup_average_cov_theta', ...
			'_all_beta_coup_average_corr_theta'};
		
		titles = {{['causal emergence ccs micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence mmi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence ccs micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence mmi micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence ccs micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence mmi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence ccs micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence mmi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence ccs micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence mmi micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical pairwise synchrony micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical sigma chi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical pairwise synchrony micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical sigma chi micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical pairwise synchrony micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical sigma chi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical pairwise synchrony micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal emergence practical sigma chi  micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation ccs micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation mmi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation ccs micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation mmi micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation ccs micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation mmi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation ccs micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation mmi micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical pairwise synchrony micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['downward causation practical sigma chi  micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling ccs micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling mmi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling ccs micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling mmi micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling ccs micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling mmi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling ccs micro bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling mmi bin sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi micro theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi micro sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi micro sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical pairwise synchrony micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['causal decoupling practical sigma chi  micro cos theta'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['average covariance theta'] ['npoints = ' num2str(npoints)]}, ...
			{['average correlation theta'] ['npoints = ' num2str(npoints)]}};

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

			ylabel('A');
			xlabel('beta');

			title(titles{i});
			exportgraphics(gcf, [pathout_plots network file_names{i} '_' num2str(npoints) '_' num2str(tau) '.png']);
			
		end
		%}
		
		close all;

	end 
end 

%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 8-node kuramoto oscillators, with fixed A, and varying beta

A = [2, 4, 6, 8, 10];

for p = 1:length(A);

g = A(p);

ce_ccs_temp = ce_ccs.ce_ccs(g,:);
figure;
scatter(beta_vec, ce_ccs_temp , 60, 'filled');
title(['causal emergence ccs, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergemce ccs');
xlabel('beta');

a_string = num2str(coupling_vec(g));
a_string = a_string(3:end);
exportgraphics(gcf, [pathout_plots_phiid_ce network '_ce_ccs_' a_string '.png']);

ce_mmi_temp = emergence_mmi.ce_mmi(g,:);
figure;
scatter(beta_vec, ce_mmi_temp, 60, 'filled');
title(['causal emergence mmi, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergence mmi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_phiid_ce network '_ce_mmi_' a_string '.png']);

ce_pract_pair_sync_micro_theta_temp = ce_practical.ce_pract_pair_sync_micro_theta(g,:);
figure;
scatter(beta_vec, ce_pract_pair_sync_micro_theta_temp, 60, 'filled');
title(['causal emergence practical, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergence practical pairwise synchrony');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_practical_ce network '_ce_pract_pair_sync_micro_theta_' a_string '.png']);

ce_pract_sigma_chi_micro_theta_temp = ce_practical.ce_pract_sigma_chi_micro_theta(g,:);
figure;
scatter(beta_vec, ce_pract_sigma_chi_micro_theta_temp, 60, 'filled');
title(['causal emergence practical sigma chi, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergence practical sigma chi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_pract_ce network '_ce_pract_sigma_chi_micro_theta_' a_string '.png']);

load([pathout_data_practical_ce network '_synchronies_ '  a_string '.mat'], 'synchronies');
%load([PATHOUT3 network '_synchronies_ '  '21313' '_' sim_index '.mat'], 'synchronies');

% metastability
sigma_met_mean = [];
for i = 1:length(beta_vec);
	sigma_met = squeeze(synchronies(i,:,:));
	sigma_met_mean(i) = mean(var(sigma_met'));
end

% chimera states
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
