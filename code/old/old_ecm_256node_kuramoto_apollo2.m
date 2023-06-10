%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?
% - add integrated information measures?

%% 256-NODE KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% This script implements practical synergy capacity for 256-node Kuramoto oscillators with different couplings and phase lags, two different macro variables 
% (variance of synchronies & global average pairwise synchrony between communities), and three different micro variables (thetas, cos(thetas), synchronies, and 
% binarized synchronies).

% major sections in this script:
% choice of parameters (time-lag, data length, and thresholds)
% create coupling matrices & noise correlation vectors
% simulate Kuramoto models - including all micro and macro variables - for all values of A and all values of beta
% calculate average covariance & average correlation of micro variables
% binarize micro and macro variables 
% calculate practical measures

clear all;
clear java;
close all;
clc;

cd '/its/home/ns508/job_templates/batch_serial'
addpath '/its/home/ns508/job_templates/batch_serial/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

% directories for generated data
pathout_data = ['/its/home/ns508/job_templates/batch_serial/output/analyses/'];

%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
all_npoints = [10000]; %, 10000];
taus = [1 10 100]; % 10 100]; %, 10, 100];

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
bin_threshold_sync = 0.9;
bin_threshold_pair_sync = 0.9;
bin_threshold_sigma_chi = 0.25;

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

%% simulate Kuramoto models - including all micro and macro variables - for all values of A and all values of beta

% sim_method() obtains variables 'phase', 'sigma_chi', and 'synchrony'
% synchronies for different values of beta and given value of A are stored and saved in 'synchronies'
% 'grand_mean_pair_sync' and 'raw_values' are derived using 'synchrony' and 'phase', respectively

% {

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	disp(q)
	
	for i = 1:size(coupling_matrices, 3);
		disp(i)
		coupling_matrix = coupling_matrices(:,:,i);

		coupling_string = num2str(coupling_vec(i));
		% only if coupling is not zero, we use letter sequence as of third letter in coupling_string
		if coupling_vec(i) ~= 0
			coupling_string = coupling_string(3:end);
		end
		
		% if coupling_string is only one number as string, add zero
		if length(coupling_string) == 1;
			coupling_string = [coupling_string, '0'];
		end
		
		for j = 1:length(beta_vec)
			beta = beta_vec(j);
			
			beta_string = num2str(beta);
			% only if beta is not zero, we use letter sequence as of third letter in beta_string
			if beta ~= 0
				beta_string = beta_string(3:end);
			end
			
			% if beta_string is only one number as string, add zero
			if length(beta_string) == 1;
				beta_string = [beta_string, '0'];
			end
				
			[phase, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);	 
			
			% rows: betas; columns: communities; 3rd dimension: time-points
			synchronies(j,:,:) = synchrony;

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
			
			% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
			save([pathout_data network '_phase_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'phase');
			save([pathout_data network '_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'sigma_chi');
			save([pathout_data network '_synchrony_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'synchrony');
			save([pathout_data network '_mean_pair_sync_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'grand_mean_pair_sync');
			save([pathout network '_raw_signal_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'raw_signal');

		end 
		a_string = num2str(coupling_vec(i));
		save([pathout_data network '_synchronies_'  a_string(3:end) '_' num2str(npoints) '.mat'], ...
			'synchronies');
		
		clear synchrony;
		clear synchronies;
		
	end 
end 

%}

%% calculate average covariance & average correlation of micro variables

% {

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	
	for i = 1:size(coupling_matrices, 3);
		disp(i)
		coupling_string = num2str(coupling_vec(i));
		% only if coupling is not zero, we use letter sequence as of third letter in coupling_string
		if coupling_vec(i) ~= 0
			coupling_string = coupling_string(3:end);
		end
		
		% if coupling_string is only one number as string, add zero
		if length(coupling_string) == 1;
			coupling_string = [coupling_string, '0'];
		end
		
		for j = 1:length(beta_vec)
			beta = beta_vec(j);
			
			beta_string = num2str(beta);
			% only if beta is not zero, we use letter sequence as of third letter in beta_string
			if beta ~= 0
				beta_string = beta_string(3:end);
			end
			
			% if beta_string is only one number as string, add zero
			if length(beta_string) == 1;
				beta_string = [beta_string, '0'];
			end
			
			% load simulated model with given A and beta:
			% micro variables - phases, raw signal, synchronies
			load([pathout_data network '_phase_'  coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'phase');
			load([pathout_data network '_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'synchrony');
			load([pathout_data network '_raw_signal_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'raw_signal');
			
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
	end
end

% save covariances and correlations for all values of A and all values of beta; saved filenames consist of
% network name + '_all_mean_' + 'corr_' or 'cov_' + micro variable name + number of datapoints
save([pathout_data network '_all_mean_corr_phase_' num2str(npoints) '.mat'], ...
	'all_mean_corr_phase');
save([pathout_data network '_all_mean_cov_phase_' num2str(npoints) '.mat'], ...
	'all_mean_cov_phase');

save([pathout_data network '_all_mean_corr_raw_signal_' num2str(npoints) '.mat'], ...
	'all_mean_corr_raw_signal');
save([pathout_data network '_all_mean_cov_raw_signal_' num2str(npoints) '.mat'], ...
	'all_mean_cov_raw_signal');

save([pathout_data network '_all_mean_corr_synchronies_' num2str(npoints) '.mat'], ...
	'all_mean_corr_synchronies');
save([pathout_data network '_all_mean_cov_synchronies_' num2str(npoints) '.mat'], ...
	'all_mean_cov_synchronies');

%}

%% binarize micro and macro variables

% {

set(0,'DefaultFigureVisible','off');
nbins = 100;

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	
	rng(1);
	for i = 1:size(coupling_matrices, 3);
		disp(i)
		coupling_string = num2str(coupling_vec(i));
		% only if coupling is not zero, we use letter sequence as of third letter in coupling_string
		if coupling_vec(i) ~= 0
			coupling_string = coupling_string(3:end);
		end
		
		% if coupling_string is only one number as string, add zero
		if length(coupling_string) == 1;
			coupling_string = [coupling_string, '0'];
		end
		
		for j = 1:length(beta_vec)
			
			beta_string = num2str(beta_vec(j));
			% only if beta is not zero, we use letter sequence as of third letter in beta_string
			if beta_vec(j) ~= 0
				beta_string = beta_string(3:end);
			end
			
			% if beta_string is only one number as string, add zero
			if length(beta_string) == 1;
				beta_string = [beta_string, '0'];
			end
			
			% load simulated model with given A and beta:
			% micro variables - phases, raw signal, synchronies;
			% macro variables for practical measures for causal emergence - 
			% variance of synchronies & global average pairwise synchrony
			% between communities
			load([pathout_data network '_phase_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'phase');
			load([pathout_data network '_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
				'sigma_chi');
			load([pathout_data network '_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'synchrony');
			load([pathout_data network '_raw_signal_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'raw_signal');
			load([pathout_data network '_mean_pair_sync_' coupling_string '_' beta_string '_' num2str(npoints) '.mat'], ...
				'grand_mean_pair_sync');
			
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
			save([pathout_data network '_bin_phase_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_phase');
			save([pathout_data network '_bin_raw_signal_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_raw_signal');
			save([pathout_data network '_bin_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
				'bin_sigma_chi');
			save([pathout_data network '_bin_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_synchrony');
			save([pathout_data network '_bin_pair_sync_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
				'bin_grand_mean_pair_sync');
			
		end
	end
	
	clear bin_phase;
	clear bin_synchrony;
	clear raw_signal;
	clear bin_grand_mean_pair_sync;
	clear sigma_chi;
	
end 
%}

%% calculate practical measures

% {

for q = 1:length(all_npoints);
	npoints = all_npoints(q);
	
	for z = 1:length(taus);
		tau = taus(z);
	
		rng(1);
		for i = 1:size(coupling_matrices, 3);
			disp(i)
			
			coupling_string = num2str(coupling_vec(i));
			% only if coupling is not zero, we use letter sequence as of third letter in coupling_string
			if coupling_vec(i) ~= 0
				coupling_string = coupling_string(3:end);
			end
			
			% if coupling_string is only one number as string, add zero
			if length(coupling_string) == 1;
				coupling_string = [coupling_string, '0'];
			end
			
			for j = 1:length(beta_vec)
				beta = beta_vec(j);
				
				beta_string = num2str(beta);
				% only if beta is not zero, we use letter sequence as of third letter in beta_string
				if beta ~= 0
					beta_string = beta_string(3:end);
				end
				
				% if beta_string is only one number as string, add zero
				if length(beta_string) == 1;
					beta_string = [beta_string, '0'];
				end
				
				% load simulated model with given A and beta
				load([pathout_data network '_bin_phase_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_phase');
				load([pathout_data network '_bin_raw_signal_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_raw_signal');
				load([pathout_data network '_bin_synchrony_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_synchrony');
				load([pathout_data network '_bin_pair_sync_' coupling_string '_' beta_string  '_' num2str(npoints) '.mat'], ...
					'bin_grand_mean_pair_sync');
				load([pathout_data network '_bin_sigma_chi_' coupling_string '_' beta_string '_' num2str(npoints)  '.mat'], ...
					'bin_sigma_chi');
				
				% micro variables: phases, raw signal, synchronies, binarized synchronies
				
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

		%% storing practical measures for different micro & macro variables

		[ce_practical] = store_pract_ce_in_struct_2x4(...
			ce_pract_sigma_chi_phase, ce_pract_sigma_chi_raw_signal, ce_pract_sigma_chi_sync, ce_pract_sigma_chi_sync,...
			dc_pract_sigma_chi_phase, dc_pract_sigma_chi_raw_signal, dc_pract_sigma_chi_sync, dc_pract_sigma_chi_sync, ...
			cd_pract_sigma_chi_phase, cd_pract_sigma_chi_raw_signal, cd_pract_sigma_chi_sync, cd_pract_sigma_chi_sync,...
			ce_pract_pair_sync_phase, ce_pract_pair_sync_raw_signal, ce_pract_pair_sync_sync, ce_pract_pair_sync_sync,...
			dc_pract_pair_sync_phase, dc_pract_pair_sync_raw_signal, dc_pract_pair_sync_sync, dc_pract_pair_sync_sync,...
			cd_pract_pair_sync_phase, cd_pract_pair_sync_raw_signal, cd_pract_pair_sync_sync, cd_pract_pair_sync_sync);
		
		% saved filenames consist of
		% network name + type of causal emergence + number of datapoints + time-lag
		save([pathout_data network '_ce_pract_' num2str(npoints) '_' num2str(tau) '.mat'], ...
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
		clear ce_pract_sigma_chi_raw_signal;
		clear ce_pract_pair_sync_raw_signal;
		clear dc_pract_sigma_chi_raw_signal;
		clear dc_pract_pair_sync_raw_signal;
		clear cd_pract_sigma_chi_raw_signal;
		clear cd_pract_pair_sync_raw_signal;
		
	end
	
end 

%}
