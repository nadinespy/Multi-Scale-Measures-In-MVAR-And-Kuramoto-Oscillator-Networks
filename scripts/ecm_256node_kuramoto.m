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

% directories for generated data
pathout_data = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/'];
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/sim_time_series/'];
pathout_data_sync = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/synchronies/'];
pathout_data_bin_sync = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/binarized_synchronies/'];
pathout_data_phiid = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/phiid/'];
pathout_data_pract_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/practical_ce/'];
pathout_data_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/phiid_ce/'];
pathout_data_mean_corr = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/average_corr/'];
pathout_data_mean_cov = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/average_cov/'];

% directories for generated plots
pathout_plots = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/'];
pathout_plots_pract_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/practical_ce/'];
pathout_plots_pract_ce_pair_sync = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/practical_ce/mean_pair_sync/'];
pathout_plots_pract_ce_sigma_chi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/practical_ce/sigma_chi/'];
pathout_plots_sigma_met = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/sigma_met/'];
pathout_plots_sigma_chi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/sigma_chi/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/distributions/'];
pathout_plots_mean_corr = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto/mean_corr/'];


%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
all_npoints = [2000]; %, 10000];
taus = [1]; %, 10, 100];

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_random(), chimera_metastable_model())
sim_method = @chimera_metastable_model;

% network - options: 
% (1) '2node_mvar' for 2-node network with 100 different coupling strengths & noise correlations 
% (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations 
% (if choosing sim_index = 3);
% 
% (2) '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations 
% (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations 
% (if choosing sim_index = 3));
%
% (3) '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations 
% (if choosing sim_index = 1 or 2)
% 
% (4) '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations 
% (if choosing sim_index = 1 or 2)
%
% (5) '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas 
% (if choosing sim_index = 4) 
%
% (6) '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas 
% (if choosing sim_index = 5) 

network = '256node_kuramoto';
bin_threshold = 0.8;

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
	
coupling_vec = linspace(0.08, 0.8, 10);
beta_vec = linspace(0.04, 0.4, 10);				   % use beta values only up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
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
			if i~=j															    % ignore diagonals
				y1 = mod(ceil(j/intra_comm_size)-1,n_communities)+1;		% community number
				if x1 == y1														% same community
					p = d0/intra_comm_size;
					k = k0;
				else															     % different communities
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
%}

%% simulate Kuramoto models for all values of A and all values of beta

% {

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for i = 1:size(coupling_matrices, 3);
		disp(i)
		coupling_string = num2str(coupling_vec(i));
		coupling_string = coupling_string(3:end);
		coupling_matrix = coupling_matrices(:,:,i);
			
		for j = 1:length(beta_vec)

			beta = beta_vec(j);
			beta_string = num2str(beta);
			if j ~= 1
				beta_string = beta_string(3:end);
			end
				
			[phase, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);	 
			
			% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints  
			save([pathout_data_sim_time_series network '_phase_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'phase');
			save([pathout_data_sim_time_series network '_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'sigma_chi');
			save([pathout_data_sim_time_series network '_synchrony_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'synchrony');
			
		end 
	end 
end 
%}

%% calculate average covariance & average correlation, and store synchronies & binarized synchronies

% {

for i = 1:size(coupling_matrices, 3);
	disp(i)
	coupling_string = num2str(coupling_vec(i));
	coupling_string = coupling_string(3:end);
	coupling_matrix = coupling_matrices(:,:,i);
	
	for j = 1:length(beta_vec)
		
		beta = beta_vec(j);
		beta_string = num2str(beta);
		if j ~= 1
			beta_string = beta_string(3:end);
		end
		
		% load simulated model with given A and beta
		load([pathout_data_sim_time_series network '_phase_'  coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
			'phase');
		load([pathout_data_sim_time_series network '_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
			'sigma_chi');
		load([pathout_data_sim_time_series network '_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
			'synchrony');
		
		% rows: betas; columns: communities; 3rd dimension: time-points
		synchronies(j,:,:) = synchrony;
		
		% micro variables: phases, raw signal, synchronies, binarized synchronies
		
		% raw signal
		raw_signal = cos(phase);
		
		% binarize synchronies
		for k = 1:size(synchrony, 2);
			for l = 1:size(synchrony,1);
				if synchrony(l,k) >= bin_threshold;
					bin_synchrony(l,k) = 1;
				else bin_synchrony(l,k) = 0;
				end
			end
		end
		
		bin_synchronies(j,:,:) = bin_synchrony;
		
		% adding error to micro variables
		cov_err  = eye(N, N);
		mu = zeros(N,1);
		E = mvnrnd(mu, cov_err, npoints)';	
		raw_signal = raw_signal + E;
		phase = phase + E;
		
		cov_err  = eye(size(synchrony,1), size(synchrony,1));
		mu = zeros(size(synchrony,1),1);
		E = mvnrnd(mu, cov_err, npoints)';
		synchrony = synchrony + E;
		
		% average covariance/correlation matrix of micro variables
		cov_phase = cov(phase');
		all_mean_cov_phase(i,j) = mean(nonzeros(tril(cov_phase,-1)), 'all');
		
		corr_phase = corrcov(cov_phase);
		all_mean_corr_phase(i,j) = mean(nonzeros(tril(corr_phase,-1)), 'all');
		
		cov_raw_signal = cov(raw_signal');
		all_mean_cov_raw_signal(i,j) = mean(nonzeros(tril(cov_raw_signal,-1)), 'all');
		
		corr_raw_signal = corrcov(cov_raw_signal);
		all_mean_corr_raw_signal(i,j) = mean(nonzeros(tril(corr_raw_signal,-1)), 'all');
		
		cov_synchrony = cov(synchrony');
		all_mean_cov_synchronies(i,j) = mean(nonzeros(tril(cov_synchrony,-1)), 'all');
		
		corr_synchrony = corrcov(cov_synchrony);
		all_mean_corr_synchronies(i,j) = mean(nonzeros(tril(corr_synchrony,-1)), 'all');
				
	end
	
	% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
	a_string = num2str(coupling_vec(i));
	save([pathout_data_sync network '_synchronies_'  a_string(3:end) '_' num2str(npoints) '.mat'], ...
		'synchronies');
	save([pathout_data_bin_sync network '_bin_synchronies_' a_string(3:end) '_' num2str(npoints) '.mat'], ...
		'bin_synchronies');
	
	clear synchrony;
	clear bin_synchrony;
	clear synchronies;
	clear bin_synchronies;
	
end

% save covariances and correlations for all values of A and all values of beta; saved filenames consist of
% network name + '_all_mean_' + 'corr_' or 'cov_' + micro variable name + number of datapoints
save([pathout_data_mean_corr network '_all_mean_corr_phase_' num2str(npoints) '.mat'], ...
	'all_mean_corr_phase');
save([pathout_data_mean_cov network '_all_mean_cov_phase_' num2str(npoints) '.mat'], ...
	'all_mean_cov_phase');

save([pathout_data_mean_corr network '_all_mean_corr_raw_signal_' num2str(npoints) '.mat'], ...
	'all_mean_corr_raw_signal');
save([pathout_data_mean_cov network '_all_mean_cov_raw_signal_' num2str(npoints) '.mat'], ...
	'all_mean_cov_raw_signal');

save([pathout_data_mean_corr network '_all_mean_corr_synchronies_' num2str(npoints) '.mat'], ...
	'all_mean_corr_synchronies');
save([pathout_data_mean_cov network '_all_mean_cov_synchronies_' num2str(npoints) '.mat'], ...
	'all_mean_cov_synchronies');

%}

%% calculating information atoms & practical measures

% {
nbins = 100;

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	
	for z = 1:length(taus);
		tau = taus(z);
	
		rng(1);
		for i = 1:size(coupling_matrices, 3);
			disp(i)
			coupling_string = num2str(coupling_vec(i));
			coupling_string = coupling_string(3:end);
			coupling_matrix = coupling_matrices(:,:,i);
			
			for j = 1:length(beta_vec)
				
				beta = beta_vec(j);
				beta_string = num2str(beta);
				if j ~= 1
					beta_string = beta_string(3:end);
				end
				
				% load simulated model with given A and beta
				load([pathout_data_sim_time_series network '_phase_'  coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'phase');
				load([pathout_data_sim_time_series network '_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
					'sigma_chi');
				load([pathout_data_sim_time_series network '_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'synchrony');
				
				% micro variables: phases, raw signal, synchronies, binarized synchronies
				
				% raw signal
				raw_signal = cos(phase);
				
				% binarize synchronies
				for k = 1:size(synchrony, 2);
					for l = 1:size(synchrony,1);
						if synchrony(l,k) >= bin_threshold;
							bin_synchrony(l,k) = 1;
						else bin_synchrony(l,k) = 0;
						end
					end
				end
				
				% macro variables for practical measures for causal emergence:
				% variance of synchronies & global average pairwise synchrony
				% between communities
				
				% calculate global average pairwise synchrony between communities
				grand_mean_pair_sync = zeros(1, npoints);
				
				for h = 1:M;
					mean_pair_sync_temp = zeros(1, npoints);
					for g = 1:M;
						mean_pair_sync_temp = mean_pair_sync_temp + abs((synchrony(h,:)+synchrony(g,:))/2);
					end
					mean_pair_sync_temp = mean_pair_sync_temp/M;
					grand_mean_pair_sync = grand_mean_pair_sync + mean_pair_sync_temp;
				end
				grand_mean_pair_sync = grand_mean_pair_sync/M;
				
				% adding error to micro & macro variables
				cov_err  = eye(N, N);
				mu = zeros(N,1);
				E = mvnrnd(mu, cov_err, npoints)';	
				raw_signal = raw_signal + E;
				phase = phase + E;
		
				cov_err = eye(size(synchrony,1), size(synchrony,1));
				mu = zeros(size(synchrony,1),1);
				E = mvnrnd(mu, cov_err, npoints)';
				synchrony = synchrony + E;
				
				cov_err  = eye(size(sigma_chi,1), size(sigma_chi,1));
				mu = zeros(size(sigma_chi,1),1);
				E = mvnrnd(mu, cov_err, npoints)';
				sigma_chi = sigma_chi + E;
				grand_mean_pair_sync = grand_mean_pair_sync + E;
		
				% check distributions of a subset of parameters
				variables = {phase, raw_signal, synchrony, bin_synchrony, sigma_chi, grand_mean_pair_sync};
				titles = {'phase', 'raw signal', 'synchrony', 'binarized synchrony', 'sigma chi', 'global average pairwise synchrony'};
				filenames = {'phase', 'raw_signal', 'sync', 'bin_sync', 'sigma_chi', 'pair_sync'};
				
				if (((i == 1) && (j == 1)) || ((i == 3) && (j == 3)) || ((i == 7) && (j == 7)) || ((i == 10) && (j == 10)));
					
					for h = 1:length(variables);
						variable = variables{h};
						
						num_var = size(variable, 1);
						r = randi([1 num_var],1,3);
						
						for k = 1:length(r);
							a_string = num2str(coupling_vec(i));
							a_string = a_string(3:end);
							
							figure;
							histogram(variable(r(k),:)', nbins);
							title({['distribution of ', titles{h}]  ['var #', num2str(r(k)), ', A = ', num2str(coupling_vec(i)), ', beta = ', num2str(beta_vec(j))]});
							ylabel('frequency');
							xlabel(['values of ', titles{h}, ', var #', num2str(r(k))]);
							exportgraphics(gcf, [pathout_plots_distributions network '_distr_' filenames{h} '_' num2str(r(k)) '_' a_string '_' num2str(npoints) '.png']);
						end
						close all;
					
					end 
				end
				
				% ---------------------------------------------------------------------------------------------------------------------------------------
				% practical causal emergence
				% ---------------------------------------------------------------------------------------------------------------------------------------
				
				% practical causal emergence variable names consist of
				% measure + type of measure (here '_pract_') + macro variable name + micro variable name
				
				% practical causal emergence with PHASES as "micro" variables
				% & sigma chi and global average pairwise synchrony as macro variables;
				ce_pract_sigma_chi_phase(i,j) = EmergencePsi(phase', sigma_chi', tau);
				ce_pract_pair_sync_phase(i,j) = EmergencePsi(phase', grand_mean_pair_sync', tau);
				dc_pract_sigma_chi_phase(i,j) = EmergenceDelta(phase', sigma_chi', tau);
				dc_pract_pair_sync_phase(i,j) = EmergenceDelta(phase', grand_mean_pair_sync', tau);
				cd_pract_sigma_chi_phase(i,j) = ce_pract_sigma_chi_phase(i,j) - dc_pract_sigma_chi_phase(i,j);
				cd_pract_pair_sync_phase(i,j) = ce_pract_pair_sync_phase(i,j) - dc_pract_pair_sync_phase(i,j);
				
				% practical causal emergence with SYNCHRONIES as "micro" variables
				% & sigma chi and global average pairwise synchrony as macro variables
				ce_pract_sigma_chi_sync(i,j) = EmergencePsi(synchrony', sigma_chi', tau);
				ce_pract_pair_sync_sync(i,j) = EmergencePsi(synchrony', grand_mean_pair_sync', tau);
				dc_pract_sigma_chi_sync(i,j) = EmergenceDelta(synchrony', sigma_chi', tau);
				dc_pract_pair_sync_sync(i,j) = EmergenceDelta(synchrony', grand_mean_pair_sync', tau);
				cd_pract_sigma_chi_sync(i,j) = ce_pract_sigma_chi_sync(i,j) - dc_pract_sigma_chi_sync(i,j);
				cd_pract_pair_sync_sync(i,j) = ce_pract_pair_sync_sync(i,j) - dc_pract_pair_sync_sync(i,j);
				
				% doing the same, albeit with BINARIZED SYNCHRONIES as "micro" variables
				% & sigma chi and global average pairwise synchrony as macro variables
				ce_pract_sigma_chi_bin_sync(i,j) = EmergencePsi(bin_synchrony', sigma_chi', tau, 'discrete');
				ce_pract_pair_sync_bin_sync(i,j) = EmergencePsi(bin_synchrony', grand_mean_pair_sync', tau, 'discrete');
				dc_pract_sigma_chi_bin_sync(i,j) = EmergenceDelta(bin_synchrony', sigma_chi', tau, 'discrete');
				dc_pract_pair_sync_bin_sync(i,j) = EmergenceDelta(bin_synchrony', grand_mean_pair_sync', tau, 'discrete');
				cd_pract_sigma_chi_bin_sync(i,j) = ce_pract_sigma_chi_bin_sync(i,j) - dc_pract_sigma_chi_bin_sync(i,j);
				cd_pract_pair_sync_bin_sync(i,j) = ce_pract_pair_sync_bin_sync(i,j) - dc_pract_pair_sync_bin_sync(i,j);
				
				% practical causal emergence with RAW SIGNAL as the "true" micro variables
				% & sigma chi and global average pairwise synchrony as macro variables
				ce_pract_sigma_chi_raw_signal(i,j) = EmergencePsi(raw_signal', sigma_chi', tau);
				ce_pract_pair_sync_raw_signal(i,j) = EmergencePsi(raw_signal', grand_mean_pair_sync', tau);
				dc_pract_sigma_chi_raw_signal(i,j) = EmergenceDelta(raw_signal', sigma_chi', tau);
				dc_pract_pair_sync_raw_signal(i,j) = EmergenceDelta(raw_signal', grand_mean_pair_sync', tau);
				cd_pract_sigma_chi_raw_signal(i,j) = ce_pract_sigma_chi_raw_signal(i,j) - dc_pract_sigma_chi_raw_signal(i,j);
				cd_pract_pair_sync_raw_signal(i,j) = ce_pract_pair_sync_raw_signal(i,j) - dc_pract_pair_sync_raw_signal(i,j);
				
			end 
			
			clear synchrony;
			clear bin_synchrony;
			
		end

		%% storing practical measures for different micro & macro variables

		[ce_practical] = store_pract_ce_in_struct_2x4(...
			ce_pract_sigma_chi_phase, ce_pract_sigma_chi_raw_signal, ce_pract_sigma_chi_sync, ce_pract_sigma_chi_bin_sync,...
			dc_pract_sigma_chi_phase, dc_pract_sigma_chi_raw_signal, dc_pract_sigma_chi_sync, dc_pract_sigma_chi_bin_sync, ...
			cd_pract_sigma_chi_phase, cd_pract_sigma_chi_raw_signal, cd_pract_sigma_chi_sync, cd_pract_sigma_chi_bin_sync,...
			ce_pract_pair_sync_phase, ce_pract_pair_sync_raw_signal, ce_pract_pair_sync_sync, ce_pract_pair_sync_bin_sync,...
			dc_pract_pair_sync_phase, dc_pract_pair_sync_raw_signal, dc_pract_pair_sync_sync, dc_pract_pair_sync_bin_sync,...
			cd_pract_pair_sync_phase, cd_pract_pair_sync_raw_signal, cd_pract_pair_sync_sync, cd_pract_pair_sync_bin_sync);
		
		% saved filenames consist of
		% network name + type of causal emergence + number of datapoints + time-lag
		save([pathout_data_pract_ce network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_practical');

		clear ce_practical;
		clear ce_pract_sigma_chi_phase;
		clear ce_pract_pair_sync_phase;
		clear dc_pract_sigma_chi_phase;
		clear dc_pract_pair_sync_phase;
		clear cd_pract_sigma_chi_phase;
		clear cd_pract_pair_sync_phase;
		clear ce_pract_sigma_chi_sync;
		clear ce_pract_pair_sync_sync;
		clear dc_pract_sigma_chi_sync;
		clear dc_pract_pair_sync_sync;
		clear cd_pract_sigma_chi_sync;
		clear cd_pract_pair_sync_sync;
		clear ce_pract_sigma_chi_bin_sync;
		clear ce_pract_pair_sync_bin_sync;
		clear dc_pract_sigma_chi_bin_sync;
		clear dc_pract_pair_sync_bin_sync;
		clear cd_pract_sigma_chi_bin_sync;
		clear cd_pract_pair_sync_bin_sync;
		clear ce_pract_sigma_chi_raw_signal;
		clear ce_pract_pair_sync_raw_signal;
		clear dc_pract_sigma_chi_raw_signal;
		clear dc_pract_pair_sync_raw_signal;
		clear cd_pract_sigma_chi_raw_signal;
		clear cd_pract_pair_sync_raw_signal;
		
	end
	
end 

%}

%% plotting

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for z = 1:length(taus);
		tau = taus(z);

		load([pathout_data_pract_ce network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_practical');
		
		% axes ticks
		x_axis = {'0.04', '', '', '0.16', '', '', '0.28', '', '', '0.4'};
		y_axis = {'0.08', '', '', '0.32', '', '', '0.56', '', '', '0.8'};
		
		y_label = 'A';
		x_label = 'beta';
		
		%% synergistic capacity, downward causation, causal decoupling
		
		% {
		% heatmaps using matlab built-in function:
		
		% practical CE, DC, & CD for different micro variables, macro: mean pairwise synchrony
		variables = {ce_practical.ce_pract_pair_sync_phase, ...
			ce_practical.ce_pract_pair_sync_raw_signal, ...
			ce_practical.ce_pract_pair_sync_sync, ...
			ce_practical.ce_pract_pair_sync_bin_sync, ...			
			ce_practical.dc_pract_pair_sync_phase, ...
			ce_practical.dc_pract_pair_sync_raw_signal, ...
			ce_practical.dc_pract_pair_sync_sync, ...
			ce_practical.dc_pract_pair_sync_bin_sync, ...			
			ce_practical.cd_pract_pair_sync_phase, ...
			ce_practical.cd_pract_pair_sync_raw_signal, ...
			ce_practical.cd_pract_pair_sync_sync, ...
			ce_practical.cd_pract_pair_sync_bin_sync};
		
		file_names = {'_all_beta_coup_pract_ce_pair_sync_phase', ...
			'_all_beta_coup_pract_ce_pair_sync_raw_signal', ...
			'_all_beta_coup_pract_ce_pair_sync_sync', ...
			'_all_beta_coup_pract_ce_pair_sync_bin_sync', ...			
			'_all_beta_coup_pract_dc_pair_sync_phase', ...
			'_all_beta_coup_pract_dc_pair_sync_raw_signal', ...
			'_all_beta_coup_pract_dc_pair_sync_sync', ...
			'_all_beta_coup_pract_dc_pair_sync_bin_sync', ...		
			'_all_beta_coup_pract_cd_pair_sync_phase', ...
			'_all_beta_coup_pract_cd_pair_sync_raw_signal', ...
			'_all_beta_coup_pract_cd_pair_sync_sync', ...
			'_all_beta_coup_pract_cd_pair_sync_bin_sync'};
		
		titles = {{['practical CE, macro: mean pairwise sync, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: mean pairwise sync, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: mean pairwise sync, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: mean pairwise sync, micro: sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};

		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_pract_ce_pair_sync, tau);
		
		% practical CE, DC, & CD for different micro variables, macro: sigma chi
		variables = {ce_practical.ce_pract_sigma_chi_phase, ...
			ce_practical.ce_pract_sigma_chi_raw_signal, ...
			ce_practical.ce_pract_sigma_chi_bin_sync, ...
			ce_practical.ce_pract_sigma_chi_sync, ...			
			ce_practical.dc_pract_sigma_chi_phase, ...
			ce_practical.dc_pract_sigma_chi_raw_signal, ...
			ce_practical.dc_pract_sigma_chi_sync, ...
			ce_practical.dc_pract_sigma_chi_bin_sync, ...			
			ce_practical.cd_pract_sigma_chi_phase, ...
			ce_practical.cd_pract_sigma_chi_sync, ...
			ce_practical.cd_pract_sigma_chi_bin_sync, ...
			ce_practical.cd_pract_sigma_chi_raw_signal};
		
		file_names = {'_all_beta_coup_pract_ce_sigma_chi_phase', ...
			'_all_beta_coup_pract_ce_sigma_chi_raw_signal', ...
			'_all_beta_coup_pract_ce_sigma_chi_sync', ...
			'_all_beta_coup_pract_ce_sigma_chi_bin_sync', ...
			'_all_beta_coup_pract_dc_sigma_chi_phase', ...
			'_all_beta_coup_pract_dc_sigma_chi_raw_signal', ...
			'_all_beta_coup_pract_dc_sigma_chi_sync', ...
			'_all_beta_coup_pract_dc_sigma_chi_bin_sync', ...
			'_all_beta_coup_pract_cd_sigma_chi_phase', ...
			'_all_beta_coup_pract_cd_sigma_chi_bin_sync', ...
			'_all_beta_coup_pract_cd_sigma_chi_sync', ...
			'_all_beta_coup_pract_cd_sigma_chi_raw_signal'};
			
		titles = {{['practical CE, macro: sigma chi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: sigma chi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: sigma chi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: sigma chi, micro: sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: sync bin'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_pract_ce_sigma_chi, tau);
		%}
		
	end 
end 

% correlations of micro variables

load([pathout_data_mean_corr network '_all_mean_corr_phase_' num2str(npoints) '.mat'], ...
	'all_mean_corr_phase');
load([pathout_data_mean_corr network '_all_mean_corr_raw_signal_' num2str(npoints) '.mat'], ...
	'all_mean_corr_raw_signal');
load([pathout_data_mean_corr network '_all_mean_corr_synchronies_' num2str(npoints) '.mat'], ...
	'all_mean_corr_synchronies');

variables = {all_mean_corr_phase, ...
	all_mean_corr_raw_signal, ...
	all_mean_corr_synchronies};

file_names = {'_all_beta_coup_mean_corr_phase', ...
	'_all_beta_coup_mean_corr_raw_signal', ...
	'_all_beta_coup_mean_corr_synchronies'};

titles = {{['average correlation of phases'] ['npoints = ' num2str(npoints)]}, ...
	{['average correlation of raw signal'] ['npoints = ' num2str(npoints)]}, ...
	{['average correlation of synchronies'] ['npoints = ' num2str(npoints)]}};

plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_mean_corr);
		
%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 256-node kuramoto oscillators, with fixed A, and varying beta

% coupling_vec(20) = 0.2131
A = [2, 4, 6, 8, 10];

for p = 1:length(A);

g = A(p);

a_string = num2str(coupling_vec(g));
a_string = a_string(3:end);

ce_pract_pair_sync_phase_temp = ce_practical.ce_pract_pair_sync_phase(g,:);
figure;
scatter(beta_vec, ce_pract_pair_sync_phase_temp, 60, 'filled');
title(['causal emergence practical, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergence practical pairwise synchrony');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_practical_ce network '_ce_pract_pair_sync_phase_' a_string '.png']);

ce_pract_sigma_chi_phase_temp = ce_practical.ce_pract_sigma_chi_phase(g,:);
figure;
scatter(beta_vec, ce_pract_sigma_chi_phase_temp, 60, 'filled');
title(['causal emergence practical sigma chi, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergence practical sigma chi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_pract_ce network '_ce_pract_sigma_chi_phase_' a_string '.png']);

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
