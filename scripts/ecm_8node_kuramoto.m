%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?
% - add integrated information measures?

%% 8-NODE KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements causal emergence for 8-node Kuramoto oscillators with different couplings and phase lags, 
% two different macro variables (variance of synchronies & global average pairwise synchrony between communities), 
% and three different micro variables (phase, cos(phase), synchronies, and binarized synchronies).

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
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/']; 
pathout_data_sim_time_series = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/sim_time_series/'];
pathout_data_sync = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/synchronies/'];
pathout_data_bin_sync = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/bin_synchronies/'];
pathout_data_phiid = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/phiid/'];
pathout_data_pract_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/practical_ce/'];
pathout_data_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/phiid_ce/'];
pathout_data_mean_corr = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/mean_corr/'];
pathout_data_mean_cov = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/mean_cov/'];

% directories for generated plots
pathout_plots = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/'];
pathout_plots_phiid_double_red = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/phiid_double_red/'];
pathout_plots_pract_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/practical_ce/'];
pathout_plots_pract_ce_pair_sync = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/practical_ce/mean_pair_sync/'];
pathout_plots_pract_ce_sigma_chi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/practical_ce/sigma_chi/'];
pathout_plots_phiid_ce = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/phiid_ce/'];
pathout_plots_phiid_ce_mmi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/phiid_ce/mmi/'];
pathout_plots_phiid_ce_ccs = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/phiid_ce/ccs/'];
pathout_plots_sigma_met = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/sigma_met/'];
pathout_plots_sigma_chi = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/sigma_chi/'];
pathout_plots_distributions = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/distributions/'];
pathout_plots_mean_corr = ['/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/', ...
	'EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/mean_corr/'];


%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
all_npoints = [2000 10000]; %[2000, 10000];
taus = [1 10 100]; %[1, 10, 100];

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

network = '8node_kuramoto';
bin_threshold_sync = 0.9;
bin_threshold_pair_sync = 0.9;
bin_threshold_sigma_chi = 0.25;

% general order of filenames: network name, variable names, value of A, value of beta, number of datapoints, time-lag

%% load files (if already existent, to, e. g., only create plots)

%{

load([pathout_data network '_ce_phiid_ccs' sim_index '.mat'], 'ce_phiid_ccs');
load([pathout_data network '_ce_phiid_mmi' sim_index '.mat'], 'ce_phiid_mmi');
load([pathout_data network '_ce_pract ' sim_index '.mat'], 'ce_pract ');
load([pathout_data network '_phiid_all_beta_coup_ccs' sim_index '.mat'], 'phiid_all_beta_coup_ccs');
load([pathout_data network '_phiid_all_beta_coup_mmi' sim_index '.mat'], 'phiid_all_beta_coup_mmi');
load([pathout_data network '_all_mean_corr_X' sim_index '.mat'], 'all_mean_corr_X');
load([pathout_data network '_all_mean_cov_X' sim_index '.mat'], 'all_mean_cov_X');

synergy_capacity_mmi = ce_phiid_mmi.synergy_capacity_mmi;
dc_phiid_mmi = ce_phiid_mmi.dc_phiid_mmi;
cd_phiid_mmi = ce_phiid_mmi.cd_phiid_mmi;

synergy_capacity_ccs = ce_phiid_ccs.synergy_capacity_ccs;
dc_phiid_ccs = ce_phiid_ccs.dc_phiid_ccs;
cd_phiid_ccs = ce_phiid_ccs.cd_phiid_ccs;

synergy_capacity_pract _sigma_chi = ce_pract .synergy_capacity_pract _sigma_chi; 
dc_pract _sigma_chi = ce_pract .dc_pract _sigma_chi; 
cd_pract _sigma_chi = ce_pract .cd_pract _sigma_chi;

synergy_capacity_pract _pairwise_synchrony = ce_pract .synergy_capacity_pract _pairwise_synchrony; 
dc_pract _pairwise_synchrony = ce_pract .dc_pract _pairwise_synchrony; 
cd_pract _pairwise_synchrony = ce_pract .cd_pract _pairwise_synchrony;

%}

%% create coupling matrices & noise correlation vectors

% {

intra_comm_size = 4;							 % intra-community size
n_communities = 2;							  % number of communities		
	
coupling_vec = linspace(0.08, 0.8, 10);
beta_vec = linspace(0.04, 0.4, 10);				 % use beta values only up to 0.4, as sigma met & sigma chi turn out to be zero for greater 
												    % values of beta; in these cases, sigma chi will be a non-varying zero macro variable,
												    % yielding erroneous values for emergence 	
	
d0 = intra_comm_size; 
d1 = intra_comm_size;							 % numbers of connections at different community levels
		
N = intra_comm_size*n_communities;			 % total number of oscillators: 8
M = n_communities;							 % number of lowest level communities (what's that?): 2
	
coupling_matrices = zeros(N,N,length(coupling_vec));
for o = 1:length(coupling_vec);
	A = coupling_vec(o);						   % was 0.2 (the higher A, the stronger the intra-community coupling strength)
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
				else															    % different communities
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

%% calculate average covariance & average correlation, and store synchronies

% {
for q = 1:length(all_npoints);
	npoints = all_npoints(q);

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
			
			% 		% adding error to micro variables
			% 		cov_err  = eye(N, N);
			% 		mu = zeros(N,1);
			% 		E = mvnrnd(mu, cov_err, npoints)';
			% 		raw_signal = raw_signal + E;
			% 		phase = phase + E;
			%
			% 		cov_err  = eye(size(synchrony,1), size(synchrony,1));
			% 		mu = zeros(size(synchrony,1),1);
			% 		E = mvnrnd(mu, cov_err, npoints)';
			% 		synchrony = synchrony + E;
			
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
			
			% calculate average mutual information as opposed to correlation?
			
		end
		
		% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
		a_string = num2str(coupling_vec(i));
		save([pathout_data_sync network '_synchronies_'  a_string(3:end) '_' num2str(npoints) '.mat'], ...
			'synchronies');
		
		clear synchrony;
		clear synchronies;
		
	end
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

%% binarizing variables & plotting distributions

set(0,'DefaultFigureVisible','off');
nbins = 100;

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	
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
			load([pathout_data_sim_time_series network '_phase_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'phase');
			load([pathout_data_sim_time_series network '_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
				'sigma_chi');
			load([pathout_data_sim_time_series network '_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'synchrony');
			
			% micro variables: phases, raw signal, synchronies
			
			% raw signal
			raw_signal = cos(phase);
			
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
			
			% check distributions of a subset of parameters
			variables = {phase, raw_signal, synchrony, sigma_chi, grand_mean_pair_sync};
			titles = {'phase', 'raw signal', 'synchrony', 'sigma chi', 'global average pairwise synchrony'};
			filenames = {'phase', 'raw_signal', 'sync', 'sigma_chi', 'pair_sync'};
			
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
					
						figure;
						plot(variable(r(k),1800:2000));
						title({['time-series of ', titles{h}]  ['var #', num2str(r(k)), ', A = ', num2str(coupling_vec(i)), ', beta = ', num2str(beta_vec(j))]});
						ylabel('value');
						xlabel('time');
						exportgraphics(gcf, [pathout_plots_distributions network '_time_series_' filenames{h} '_' num2str(r(k)) '_' a_string '_' num2str(npoints) '.png']);
					
					end
					close all;
					
				end
			end
			
			% binarize all variables
			
			% binarize phase & raw signal
			for k = 1:size(phase, 2);
				for l = 1:size(phase,1);
					if phase(l,k) >= 0;
						bin_phase(l,k) = 1;
					else bin_phase(l,k) = 0;
					end
					
					if raw_signal(l,k) >= 0;
						bin_raw_signal(l,k) = 1;
					else bin_raw_signal(l,k) = 0;
					end
				end
			end
			
			% binarize synchronies
			for k = 1:size(synchrony, 2);
				for l = 1:size(synchrony,1);
					if synchrony(l,k) >= bin_threshold_sync;
						bin_synchrony(l,k) = 1;
					else bin_synchrony(l,k) = 0;
					end
				end
			end
			
			% binarize sigma chi & global average pairwise synchrony
			for k = 1:size(sigma_chi, 2);
				if sigma_chi(k) >= bin_threshold_sigma_chi;
					bin_sigma_chi(k) = 1;
				else bin_sigma_chi(k) = 0;
				end
				if grand_mean_pair_sync(k) >= bin_threshold_pair_sync;
					bin_grand_mean_pair_sync(k) = 1;
				else bin_grand_mean_pair_sync(k) = 0;
				end
			end
			
			% load simulated model with given A and beta
			save([pathout_data_sim_time_series network '_bin_phase_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_phase');
			save([pathout_data_sim_time_series network '_bin_raw_signal_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_raw_signal');
			save([pathout_data_sim_time_series network '_bin_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
				'bin_sigma_chi');
			save([pathout_data_sim_time_series network '_bin_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_synchrony');
			save([pathout_data_sim_time_series network '_bin_pair_sync_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_grand_mean_pair_sync');
			
		end
	end
	
	clear bin_phase;
	clear bin_synchrony;
	clear raw_signal;
	clear bin_grand_mean_pair_sync;
	clear sigma_chi;
	
end 

				
%% calculating information atoms & practical measures

% {

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
				load([pathout_data_sim_time_series network '_bin_phase_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_phase');
				load([pathout_data_sim_time_series network '_bin_raw_signal_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_raw_signal');
				load([pathout_data_sim_time_series network '_bin_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_synchrony');
				load([pathout_data_sim_time_series network '_bin_pair_sync_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_grand_mean_pair_sync');
				load([pathout_data_sim_time_series network '_bin_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
					'bin_sigma_chi');
				
% 				% adding error to micro & macro variables
% 				cov_err  = eye(N, N);
% 				mu = zeros(N,1);
% 				E = mvnrnd(mu, cov_err, npoints)';	
% 				raw_signal = raw_signal + E;
% 				phase = phase + E;
% 		
% 				cov_err = eye(size(synchrony,1), size(synchrony,1));
% 				mu = zeros(size(synchrony,1),1);
% 				E = mvnrnd(mu, cov_err, npoints)';
% 				synchrony = synchrony + E;
% 				
% 				cov_err  = eye(size(sigma_chi,1), size(sigma_chi,1));
% 				mu = zeros(size(sigma_chi,1),1);
% 				E = mvnrnd(mu, cov_err, npoints)';
% 				sigma_chi = sigma_chi + E;
% 				grand_mean_pair_sync = grand_mean_pair_sync + E;
				
% 				% take summation of phases and synchrony as (meaningless) macro variables
% 				macro_var_sum1 = zeros(1, npoints);
% 				for k = 1:(size(phase,1));
% 					macro_var_sum1 = macro_var_sum1 + phase(k,:);
% 				end
% 				
% 				macro_var_sum2 = zeros(1, npoints);
% 				for k = 1:(size(synchrony,1));
% 					macro_var_sum2 = macro_var_sum2 + synchrony(k,:);
% 				end
% 				
% 				sigma_chi = macro_var_sum1;
% 				grand_mean_pair_sync = macro_var_sum2;

				
				% ---------------------------------------------------------------------------------------------------------------------------------------
				% integrated information decomposition (PhiID)
				% ---------------------------------------------------------------------------------------------------------------------------------------
				
				% PhiID variable names consist of
				% 'phiid_' + parameters involved (here 'all_beta_coup_') + redundancy function + micro variable name
				
				% do PhiID with PHASES as micro variables
				try
					phiid_all_beta_coup_mmi_phase(:,i,j) = struct2array(PhiIDFullDiscrete(bin_phase, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_phase(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_phase(:,i,j) = struct2array(PhiIDFullDiscrete(bin_phase, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_phase(:,i,j) = NaN;
				end
				
				% do PhiID with RAW SIGNAL as micro variables
				try
					phiid_all_beta_coup_mmi_raw_signal(:,i,j) = struct2array(PhiIDFullDiscrete(bin_raw_signal, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_raw_signal(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_raw_signal(:,i,j) = struct2array(PhiIDFullDiscrete(bin_raw_signal, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_raw_signal(:,i,j) = NaN;
				end
				
				% do PhiID with SYNCHRONIES as micro variables
				try
					phiid_all_beta_coup_mmi_sync(:,i,j) = struct2array(PhiIDFullDiscrete(bin_synchrony, tau, 'MMI'))';
				catch
					phiid_all_beta_coup_mmi_sync(:,i,j) = NaN;
				end
				
				try
					phiid_all_beta_coup_ccs_sync(:,i,j) = struct2array(PhiIDFullDiscrete(bin_synchrony, tau, 'ccs'))';
				catch
					phiid_all_beta_coup_ccs_sync(:,i,j) = NaN;
				end

				
				% ---------------------------------------------------------------------------------------------------------------------------------------
				% practical causal emergence
				% ---------------------------------------------------------------------------------------------------------------------------------------
				
				% practical causal emergence variable names consist of
				% measure + type of measure (here '_pract_') + macro variable name + micro variable name
				
				% practical causal emergence with PHASES as "micro" variables
				% & sigma chi and global average pairwise synchrony as macro variables;
				ce_pract_sigma_chi_phase(i,j) = EmergencePsi(bin_phase', bin_sigma_chi', tau, 'discrete');
				ce_pract_pair_sync_phase(i,j) = EmergencePsi(bin_phase', bin_grand_mean_pair_sync', tau, 'discrete');
				dc_pract_sigma_chi_phase(i,j) = EmergenceDelta(bin_phase', bin_sigma_chi', tau, 'discrete');
				dc_pract_pair_sync_phase(i,j) = EmergenceDelta(bin_phase', bin_grand_mean_pair_sync', tau, 'discrete');
				cd_pract_sigma_chi_phase(i,j) = ce_pract_sigma_chi_phase(i,j) - dc_pract_sigma_chi_phase(i,j);
				cd_pract_pair_sync_phase(i,j) = ce_pract_pair_sync_phase(i,j) - dc_pract_pair_sync_phase(i,j);
				
				% practical causal emergence with SYNCHRONIES as "micro" variables
				% & sigma chi and global average pairwise synchrony as macro variables
				ce_pract_sigma_chi_sync(i,j) = EmergencePsi(bin_synchrony', bin_sigma_chi', tau, 'discrete');
				ce_pract_pair_sync_sync(i,j) = EmergencePsi(bin_synchrony', bin_grand_mean_pair_sync', tau, 'discrete');
				dc_pract_sigma_chi_sync(i,j) = EmergenceDelta(bin_synchrony', bin_sigma_chi', tau, 'discrete');
				dc_pract_pair_sync_sync(i,j) = EmergenceDelta(bin_synchrony', bin_grand_mean_pair_sync', tau, 'discrete');
				cd_pract_sigma_chi_sync(i,j) = ce_pract_sigma_chi_sync(i,j) - dc_pract_sigma_chi_sync(i,j);
				cd_pract_pair_sync_sync(i,j) = ce_pract_pair_sync_sync(i,j) - dc_pract_pair_sync_sync(i,j);
				
				% practical causal emergence with RAW SIGNAL as the "true" micro variables
				% & sigma chi and global average pairwise synchrony as macro variables
				ce_pract_sigma_chi_raw_signal(i,j) = EmergencePsi(bin_raw_signal', bin_sigma_chi', tau, 'discrete');
				ce_pract_pair_sync_raw_signal(i,j) = EmergencePsi(bin_raw_signal', bin_grand_mean_pair_sync', tau, 'discrete');
				dc_pract_sigma_chi_raw_signal(i,j) = EmergenceDelta(bin_raw_signal', bin_sigma_chi', tau, 'discrete');
				dc_pract_pair_sync_raw_signal(i,j) = EmergenceDelta(bin_raw_signal', bin_grand_mean_pair_sync', tau, 'discrete');
				cd_pract_sigma_chi_raw_signal(i,j) = ce_pract_sigma_chi_raw_signal(i,j) - dc_pract_sigma_chi_raw_signal(i,j);
				cd_pract_pair_sync_raw_signal(i,j) = ce_pract_pair_sync_raw_signal(i,j) - dc_pract_pair_sync_raw_signal(i,j);
				
			end 
			
			clear synchrony;
			clear raw_signal;
			clear phase;
			clear grand_mean_pair_sync;
			clear sigma_chi;
			
		end
		
		%% storing practical measures for different micro & macro variables, and information atoms in struct files
		
		[phiid_all_beta_coup_mmi_phase, phiid_all_beta_coup_ccs_phase] =  ...
			store_atoms_in_struct(phiid_all_beta_coup_mmi_phase, phiid_all_beta_coup_ccs_phase);
		[phiid_all_beta_coup_mmi_raw_signal, phiid_all_beta_coup_ccs_raw_signal] =  ...
			store_atoms_in_struct(phiid_all_beta_coup_mmi_raw_signal, phiid_all_beta_coup_ccs_raw_signal);
		[phiid_all_beta_coup_mmi_sync, phiid_all_beta_coup_ccs_sync] =  ...
			store_atoms_in_struct(phiid_all_beta_coup_mmi_sync, phiid_all_beta_coup_ccs_sync);
		
		% saved filenames consist of
		% network name + '_phiid_' + parameters involved (here 'all_beta_coup_') + redundancy function +
		% micro variable name + number of datapoints + time-lag
		save([pathout_data_phiid network '_phiid_all_beta_coup_ccs_phase_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_ccs_phase');
		save([pathout_data_phiid network '_phiid_all_beta_coup_mmi_phase_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_mmi_phase');
		
		save([pathout_data_phiid network '_phiid_all_beta_coup_ccs_raw_signal_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_ccs_raw_signal');
		save([pathout_data_phiid network '_phiid_all_beta_coup_mmi_raw_signal_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_mmi_raw_signal');
		
		save([pathout_data_phiid network '_phiid_all_beta_coup_ccs_sync_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_ccs_sync');
		save([pathout_data_phiid network '_phiid_all_beta_coup_mmi_sync_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_mmi_sync');
		
		%% calculate phiid-based synergistic/emergent capacity, downward causation, causal decoupling
		
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
		
		% causal emergence with phase as micro variables
		[ce_phiid_mmi_phase, ...
			dc_phiid_mmi_phase, ...
			cd_phiid_mmi_phase, ...
			ce_phiid_ccs_phase, ...
			dc_phiid_ccs_phase, ...
			cd_phiid_ccs_phase] = phiid_ce(phiid_all_beta_coup_mmi_phase, phiid_all_beta_coup_ccs_phase);
		
		% causal emergence with raw signal as micro variables
		[ce_phiid_mmi_raw_signal, ...
			dc_phiid_mmi_raw_signal, ...
			cd_phiid_mmi_raw_signal, ...
			ce_phiid_ccs_raw_signal, ...
			dc_phiid_ccs_raw_signal, ...
			cd_phiid_ccs_raw_signal] = phiid_ce(phiid_all_beta_coup_mmi_raw_signal, phiid_all_beta_coup_ccs_raw_signal);
		
		% causal emergence with synchronies as micro variables
		[ce_phiid_mmi_sync, ...
			dc_phiid_mmi_sync, ...
			cd_phiid_mmi_sync, ...
			ce_phiid_ccs_sync, ...
			dc_phiid_ccs_sync, ...
			cd_phiid_ccs_sync] = phiid_ce(phiid_all_beta_coup_mmi_sync, phiid_all_beta_coup_ccs_sync);
		
		% save variables in a struct
		
		[ce_phiid_ccs, ce_phiid_mmi, ce_practical] = store_ce_in_struct_2x4(...
			ce_phiid_mmi_phase, ce_phiid_mmi_raw_signal, ce_phiid_mmi_sync, ce_phiid_mmi_sync,...
			dc_phiid_mmi_phase, dc_phiid_mmi_raw_signal, dc_phiid_mmi_sync, dc_phiid_mmi_sync, ...
			cd_phiid_mmi_phase, cd_phiid_mmi_raw_signal, cd_phiid_mmi_sync, cd_phiid_mmi_sync,...
			ce_phiid_ccs_phase, ce_phiid_ccs_raw_signal, ce_phiid_ccs_sync, ce_phiid_ccs_sync,...
			dc_phiid_ccs_phase, dc_phiid_ccs_raw_signal, dc_phiid_ccs_sync, dc_phiid_ccs_sync, ...
			cd_phiid_ccs_phase, cd_phiid_ccs_raw_signal, cd_phiid_ccs_sync, cd_phiid_ccs_sync,...
			ce_pract_sigma_chi_phase, ce_pract_sigma_chi_raw_signal, ce_pract_sigma_chi_sync, ce_pract_sigma_chi_sync,...
			dc_pract_sigma_chi_phase, dc_pract_sigma_chi_raw_signal, dc_pract_sigma_chi_sync, dc_pract_sigma_chi_sync, ...
			cd_pract_sigma_chi_phase, cd_pract_sigma_chi_raw_signal, cd_pract_sigma_chi_sync, cd_pract_sigma_chi_sync,...
			ce_pract_pair_sync_phase, ce_pract_pair_sync_raw_signal, ce_pract_pair_sync_sync, ce_pract_pair_sync_sync,...
			dc_pract_pair_sync_phase, dc_pract_pair_sync_raw_signal, dc_pract_pair_sync_sync, dc_pract_pair_sync_sync,...
			cd_pract_pair_sync_phase, cd_pract_pair_sync_raw_signal, cd_pract_pair_sync_sync, cd_pract_pair_sync_sync);
		
		% saved filenames consist of
		% network name + type of causal emergence + number of datapoints + time-lag
		save([pathout_data_phiid_ce network '_ce_phiid_ccs_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_phiid_ccs');
		save([pathout_data_phiid_ce network '_ce_phiid_mmi_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_phiid_mmi');
		save([pathout_data_pract_ce network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_practical');

		clear ce_phiid_ccs;
		clear ce_phiid_mmi;
		clear ce_practical;
		clear phiid_all_beta_coup_ccs_phase;
		clear phiid_all_beta_coup_mmi_phase;
		clear phiid_all_beta_coup_ccs_raw_signal;
		clear phiid_all_beta_coup_mmi_raw_signal;
		clear phiid_all_beta_coup_ccs_sync;
		clear phiid_all_beta_coup_mmi_sync;
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

set(0,'DefaultFigureVisible','off')

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	
	for z = 1:length(taus);
		tau = taus(z);
		
		load([pathout_data_phiid network '_phiid_all_beta_coup_ccs_phase_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_ccs_phase');
		load([pathout_data_phiid network '_phiid_all_beta_coup_mmi_phase_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_mmi_phase');
		load([pathout_data_phiid network '_phiid_all_beta_coup_ccs_raw_signal_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_ccs_raw_signal');
		load([pathout_data_phiid network '_phiid_all_beta_coup_mmi_raw_signal_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_mmi_raw_signal');
		load([pathout_data_phiid network '_phiid_all_beta_coup_ccs_sync_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_ccs_sync');
		load([pathout_data_phiid network '_phiid_all_beta_coup_mmi_sync_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'phiid_all_beta_coup_mmi_sync');
		
		load([pathout_data_phiid_ce network '_ce_phiid_ccs_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_phiid_ccs');
		load([pathout_data_phiid_ce network '_ce_phiid_mmi_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_phiid_mmi');
		load([pathout_data_pract_ce network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], ...
			'ce_practical');
		
		% axes ticks and labels
		x_axis = {'0.04', '', '', '0.16', '', '', '0.28', '', '', '0.4'};
		y_axis = {'0.08', '', '', '0.32', '', '', '0.56', '', '', '0.8'};
		
		y_label = 'A';
		x_label = 'beta';


		%% double-redundancy
		
		% {
		
		double_red = {phiid_all_beta_coup_ccs_phase.rtr, ...
			phiid_all_beta_coup_mmi_phase.rtr, ...
			phiid_all_beta_coup_ccs_raw_signal.rtr, ...
			phiid_all_beta_coup_mmi_raw_signal.rtr, ...
			phiid_all_beta_coup_ccs_sync.rtr, ...
			phiid_all_beta_coup_mmi_sync.rtr};
		
		file_names = {'_all_beta_coup_ccs_rtr_phase', ...
			'_all_beta_coup_mmi_rtr_phase', ...
			'_all_beta_coup_ccs_rtr_raw_signal', ...
			'_all_beta_coup_mmi_rtr_raw_signal', ...
			'_all_beta_coup_ccs_rtr_sync', ...
			'_all_beta_coup_mmi_rtr_sync'};
		
		titles = {{['double-redundancy ccs, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy ccs, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy ccs, micro: synchronies'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['double-redundancy mmi, micro: synchronies'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		plot_heatmap(double_red, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_phiid_double_red, tau);

		%}
		
		%% synergistic capacity, downward causation, causal decoupling
		
		% {
		
		% PhiID-based CE, DC, & CD for different micro variables, using CCS
		variables = {ce_phiid_ccs.ce_phiid_ccs_phase, ...
			ce_phiid_ccs.ce_phiid_ccs_raw_signal, ...
			ce_phiid_ccs.ce_phiid_ccs_sync, ...			
			ce_phiid_ccs.dc_phiid_ccs_phase, ...
			ce_phiid_ccs.dc_phiid_ccs_raw_signal, ...
			ce_phiid_ccs.dc_phiid_ccs_sync, ...		
			ce_phiid_ccs.cd_phiid_ccs_phase, ...
			ce_phiid_ccs.cd_phiid_ccs_raw_signal, ...
			ce_phiid_ccs.cd_phiid_ccs_sync};		
			
		file_names = {'_all_beta_coup_ccs_ce_phase', ...
			'_all_beta_coup_ccs_ce_raw_signal', ...
			'_all_beta_coup_ccs_ce_sync', ...			
			'_all_beta_coup_ccs_dc_phase', ...
			'_all_beta_coup_ccs_dc_raw_signal', ...
			'_all_beta_coup_ccs_dc_sync', ...		
			'_all_beta_coup_ccs_cd_phase', ...
			'_all_beta_coup_ccs_cd_raw_signal', ...
			'_all_beta_coup_ccs_cd_sync'};
		
		titles = {{['PhiID-based CE ccs, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CE ccs, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CE ccs, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based DC ccs, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based DC ccs, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based DC ccs, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CD ccs, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CD ccs, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CD ccs, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_phiid_ce_ccs, tau);
		
		% PhiID-based CE, DC, & CD for different micro variables, using MMI
		variables = {ce_phiid_mmi.ce_phiid_mmi_phase, ...
			ce_phiid_mmi.ce_phiid_mmi_raw_signal, ...
			ce_phiid_mmi.ce_phiid_mmi_sync, ...	
			ce_phiid_mmi.dc_phiid_mmi_phase, ...
			ce_phiid_mmi.dc_phiid_mmi_raw_signal, ...
			ce_phiid_mmi.dc_phiid_mmi_sync, ...			
			ce_phiid_mmi.cd_phiid_mmi_phase, ...
			ce_phiid_mmi.cd_phiid_mmi_raw_signal, ...
			ce_phiid_mmi.cd_phiid_mmi_sync};
		
		file_names = {'_all_beta_coup_mmi_ce_phase', ...
			'_all_beta_coup_mmi_ce_raw_signal', ...
			'_all_beta_coup_mmi_ce_sync', ...		
			'_all_beta_coup_mmi_dc_phase', ...
			'_all_beta_coup_mmi_dc_raw_signal', ...
			'_all_beta_coup_mmi_dc_sync', ...			
			'_all_beta_coup_mmi_cd_phase', ...
			'_all_beta_coup_mmi_cd_raw_signal', ...
			'_all_beta_coup_mmi_cd_sync'};
		
		titles = {{['PhiID-based CE mmi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CE mmi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CE mmi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based DC mmi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based DC mmi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based DC mmi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CD mmi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CD mmi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['PhiID-based CD mmi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_phiid_ce_mmi, tau);
		
		% practical CE, DC, & CD for different micro variables, macro: mean pairwise synchrony
		variables = {ce_practical.ce_pract_pair_sync_phase, ...
			ce_practical.ce_pract_pair_sync_raw_signal, ...
			ce_practical.ce_pract_pair_sync_sync, ...		
			ce_practical.dc_pract_pair_sync_phase, ...
			ce_practical.dc_pract_pair_sync_raw_signal, ...
			ce_practical.dc_pract_pair_sync_sync, ...		
			ce_practical.cd_pract_pair_sync_phase, ...
			ce_practical.cd_pract_pair_sync_raw_signal, ...
			ce_practical.cd_pract_pair_sync_sync};
		
		file_names = {'_all_beta_coup_pract_ce_pair_sync_phase', ...
			'_all_beta_coup_pract_ce_pair_sync_raw_signal', ...
			'_all_beta_coup_pract_ce_pair_sync_sync', ...			
			'_all_beta_coup_pract_dc_pair_sync_phase', ...
			'_all_beta_coup_pract_dc_pair_sync_raw_signal', ...
			'_all_beta_coup_pract_dc_pair_sync_sync', ...		
			'_all_beta_coup_pract_cd_pair_sync_phase', ...
			'_all_beta_coup_pract_cd_pair_sync_raw_signal', ...
			'_all_beta_coup_pract_cd_pair_sync_sync'};
		
		titles = {{['practical CE, macro: mean pairwise sync, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: mean pairwise sync, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: mean pairwise sync, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: mean pairwise sync, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: mean pairwise sync, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_pract_ce_pair_sync, tau);
		
		% practical CE, DC, & CD for different micro variables, macro: sigma chi
		variables = {ce_practical.ce_pract_sigma_chi_phase, ...
			ce_practical.ce_pract_sigma_chi_raw_signal, ...
			ce_practical.ce_pract_sigma_chi_sync, ...			
			ce_practical.dc_pract_sigma_chi_phase, ...
			ce_practical.dc_pract_sigma_chi_raw_signal, ...
			ce_practical.dc_pract_sigma_chi_sync, ...		
			ce_practical.cd_pract_sigma_chi_phase, ...
			ce_practical.cd_pract_sigma_chi_sync, ...
			ce_practical.cd_pract_sigma_chi_raw_signal};
		
		file_names = {'_all_beta_coup_pract_ce_sigma_chi_phase', ...
			'_all_beta_coup_pract_ce_sigma_chi_raw_signal', ...
			'_all_beta_coup_pract_ce_sigma_chi_sync', ...
			'_all_beta_coup_pract_dc_sigma_chi_phase', ...
			'_all_beta_coup_pract_dc_sigma_chi_raw_signal', ...
			'_all_beta_coup_pract_dc_sigma_chi_sync', ...
			'_all_beta_coup_pract_cd_sigma_chi_phase', ...
			'_all_beta_coup_pract_cd_sigma_chi_sync', ...
			'_all_beta_coup_pract_cd_sigma_chi_raw_signal'};
		
		titles = {{['practical CE, macro: sigma chi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: sigma chi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CE, macro: sigma chi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical DC, macro: sigma chi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: phase'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: raw signal'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}, ...
			{['practical CD, macro: sigma chi, micro: sync'] ['npoints = ' num2str(npoints) ', tau = ' num2str(tau)]}};
		
		plot_heatmap(variables, file_names, titles, x_axis, y_axis, x_label, y_label, network, npoints, pathout_plots_pract_ce_sigma_chi, tau);

		%}
		
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

end 


%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 8-node kuramoto oscillators, with fixed A, and varying beta

A = [2, 4, 6, 8, 10];

for p = 1:length(A);

g = A(p);

ce_phiid_ccs_temp = ce_phiid_ccs.ce_phiid_ccs(g,:);
figure;
scatter(beta_vec, ce_phiid_ccs_temp , 60, 'filled');
title(['causal emergence ccs, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergemce ccs');
xlabel('beta');

a_string = num2str(coupling_vec(g));
a_string = a_string(3:end);
exportgraphics(gcf, [pathout_plots_phiid_ce network '_ce_phiid_ccs_' a_string '.png']);

ce_phiid_mmi_temp = emergence_phiid_mmi.ce_phiid_mmi(g,:);
figure;
scatter(beta_vec, ce_phiid_mmi_temp, 60, 'filled');
title(['causal emergence mmi, A = ' num2str(coupling_vec(g))]);
ylabel('causal emergence mmi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots_phiid_ce network '_ce_phiid_mmi_' a_string '.png']);

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
