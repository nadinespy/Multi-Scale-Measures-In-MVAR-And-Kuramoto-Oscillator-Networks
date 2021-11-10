%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?

% add macro-variables
% add integrated information measures?

%% 
% This script implements synergy capacity for 2-node and 8-node MVAR models with differing connection strengths and noise correlations

clear all;
clear java;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

% 2-node mvar networks with different couplings
PATHOUT1a = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/2node_mvar/';
PATHOUT1b = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/2node_mvar/';

% 8-node mvar networks with different architectures
PATHOUT2a = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_mvar_different_architectures/';
PATHOUT2b = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_mvar_different_architectures/';

% 8-node mvar networks with global coupling factor
PATHOUT3a = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_mvar_global_coupling/';
PATHOUT3b = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_mvar_global_coupling/';

% 8-node mvar ER networks
PATHOUT4a = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_mvar_erdoes_renyi/';
PATHOUT4b = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_mvar_erdoes_renyi/';

% 8-node kuramoto oscillators
PATHOUT5a = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/';
PATHOUT5b = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto';

% 256-node kuramoto oscillators
PATHOUT6a = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/256node_kuramoto/';
PATHOUT6b = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/256node_kuramoto';



%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
npoints = 2000;
tau = 1;

% ------------------------------------------------------------------------------------------------------
% MODEL
% ------------------------------------------------------------------------------------------------------

% MULTIVARIATE AUTOREGRESSIVE TIME-SERIES (MVAR)

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_random(), chimera_metastable_model())
sim_method = @chimera_metastable_model;

% save plots & matrices according to simulation method (options: '1' (for statdata_coup_errors1()), '2' (for statdata_coup_errors2()), '3' (for statdata_random(), '4' (for metastable_chimera_model())
sim_index = '4';

% network (options: '2node' for 2-node network with 100 different coupling strengths & noise correlations (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations (if choosing sim_index = 3);
% '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations (if choosing sim_index = 3));
% '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations (if choosing sim_index = 1 or 2)
% '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations (if choosing sim_index = 1 or 2)
% '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 4) 
% '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 5) 

network = '256node_kuramoto';

% choose according to which model is being investigated (see paths above)
pathout_plots = PATHOUT6a;
pathout_data = PATHOUT6b;

% KURAMOTO OSCILLATORS


% ------------------------------------------------------------------------------------------------------
% MEASURE
% ------------------------------------------------------------------------------------------------------


% ------------------------------------------------------------------------------------------------------
% MACRO VARIABLE
% ------------------------------------------------------------------------------------------------------



%% load files (if already existent, to, e. g., only create plots)

%{

if strcmp(network, '256node_kuramoto') == true;
	load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');
	synergy_capacity_practical = emergence_practical.synergy_capacity_practical; 
	downward_causation_practical = emergence_practical.downward_causation_practical; 
	causal_decoupling_practical = emergence_practical.causal_decoupling_practical;

else
	load([pathout_data network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
	load([pathout_data network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
	load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');
	load([pathout_data network '_all_atoms_err_coup_ccs' sim_index '.mat'], 'all_atoms_err_coup_ccs');
	load([pathout_data network '_all_atoms_err_coup_mmi' sim_index '.mat'], 'all_atoms_err_coup_mmi');
	%load([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
	%load([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

	synergy_capacity_mmi = emergence_mmi.synergy_capacity_mmi; 
	downward_causation_mmi = emergence_mmi.downward_causation_mmi; 
	causal_decoupling_mmi = emergence_mmi.causal_decoupling_mmi; 

	synergy_capacity_ccs = emergence_ccs.synergy_capacity_ccs; 
	downward_causation_ccs = emergence_ccs.downward_causation_ccs; 
	causal_decoupling_ccs = emergence_ccs.causal_decoupling_ccs; 

	synergy_capacity_practical = emergence_practical.synergy_capacity_practical; 
	downward_causation_practical = emergence_practical.downward_causation_practical; 
	causal_decoupling_practical = emergence_practical.causal_decoupling_practical;

end

%}

%% create coupling matrices (& noise correlations) for the different models

% {
if strcmp(network, '2node_mvar') == true; 
	
	if sim_index == '3';
		coupling_vec = linspace(0.0, 0.0, 100);
		error_vec = linspace(0.0, 0.0, 100); 
	else coupling_vec = linspace(0.01,0.45, 100);
		error_vec = linspace(0.01, 0.9, 100); 
	end 

	coupling_matrices = []; 
	for i = 1:length(coupling_vec);
		coupling_matrices(:,:,i) = coupling_vec(i)*ones(2);
	end 
	
elseif strcmp(network, '8node_mvar_different_architectures') == true; 
	if sim_index == '3'

		coupling_vec = linspace(0.0, 0.0, 100);
		error_vec = linspace(0.0, 0.0, 100);
		
		coupling_matrices = []; 
		for i = 1:length(coupling_vec);
			coupling_matrices(:,:,i) = coupling_vec(i)*ones(8);
		end 
	
	else load('all_nets.mat');									% phi-optimal binary network, phi-optimal weighted network, small world, fully connected, bidirectional ring, unidirectional ring
		
		net_names = fields(all_nets);
		
		error_vec = linspace(0.01, 0.9, 7);
		
		coupling_matrices = [];
		all_sr = [];
		for i = 1:size(net_names,1);
			net = all_nets.(net_names{i});			% normalise network so that spectral radius is close to but below 1
			sr = max(abs(eig(net)));
			all_sr.(net_names{i}) = sr;
                  net = net./(1.10*sr);
			coupling_matrices(:,:,i) = net;
		end 
		save([PATHOUT1 network '_spectral_radius' sim_index '.mat'], 'all_sr');
		
	end

elseif strcmp(network, '8node_mvar_erdoes_renyi') == true;
	
	density_vec = linspace(0.05,0.9, 100);
	error_vec = linspace(0.01, 0.9, 100);
	nb_ER_networks = 50;
	
	coupling_matrices = [];
	
	for n = 1:length(error_vec);						% to continue with PhiG: n = 7
		noise = error_vec(n);	
		
		for o = 1:length(density_vec);
			density = density_vec(o);
			
			fifty_values_per_ER_network = zeros(50, 1);
			for p = 1:nb_ER_networks

				ER_network = zeros(8,8);
				a = 0;
				b = 0.99;
				for i = 1:length(ER_network);
					for j = 1:length(ER_network);
						r = (b-a) .* rand(1) + a;
			
						if (r <= density && i ~= j)
							ER_network(i, j) = 1;
						end 
					end 
				end 
				
				% normalise network so that spectral radius is close to but below 1
				sr = max(abs(eig(ER_network)));
				ER_network_normalized = ER_network./(1.10*sr);
				coupling_matrices(:,:,o,p) = ER_network_normalized;
				
			end
		end
	end

elseif strcmp(network, '8node_mvar_global_coupling') == true;
	load('all_nets.mat');									% phi-optimal binary network, phi-optimal weighted network, small world, fully connected, bidirectional ring, unidirectional ring
	net_names = fields(all_nets);
	optimalB_net = all_nets.OptimalB;
	coupling_vec = linspace(0.05,0.9, 100);
	error_vec = linspace(0.01, 0.9, 100);
	coupling_matrices = [];
	
% optimalB_net(optimalB_net ~= 0) = 1;				% Pedro took the original weighted network as opposed to a binary one --> why should I opt for one or the other? This decision completely changes results for some measures 
	k_vec = linspace(0.01, 0.9, 100000);

	for i = 1:length(k_vec)
		sr = max(abs(eig(k_vec(i) * optimalB_net)));
		if sr >= 1
			sr = max(abs(eig(k_vec(i-1) * optimalB_net)));
			k_max = k_vec(i-1);
			break
		else k_max = k_vec(i);	
		end 
	end 
	
	for o = 1:length(coupling_vec);
			coupling_matrices(:,:,o) = k_max * coupling_vec(o) * optimalB_net;
	end 

elseif ((strcmp(network, '8node_kuramoto') == true) || (strcmp(network, '256node_kuramoto') == true));
	
	if strcmp(network, '8node_kuramoto') == true
		intra_comm_size = 4;						% intra-community size
		n_communities = 2;							% number of communities
		coupling_vec = linspace(0.05, 0.9, 100);
		error_vec = linspace(0.0, 0.8, 100);			% not an error vector, but beta (named it error_vec so that i don't need to rewrite everything below)
	else intra_comm_size = 32;
		n_communities = 8;
		coupling_vec = linspace(0.05, 0.9, 10);
		error_vec = linspace(0.0, 0.8, 10);
	end
		
			
	d0 = intra_comm_size; 
	d1 = intra_comm_size;			     % numbers of connections at different community levels
		
	N = intra_comm_size*n_communities;	% total number of oscillators: 64
	M = n_communities;				       % number of lowest level communities (what's that?): 4
	
	synchronies = zeros(length(error_vec), n_communities, npoints);
	
	coupling_matrices = zeros(N,N,length(coupling_vec));
	for o = 1:length(coupling_vec);
		A = coupling_vec(o);			% was 0.2 (the higher A, the stronger the intra-community coupling strength)
		k1 = (1-A)/2;					% inter-community coupling strength: 0.4
		k0 = 1-k1;						% intra-community coupling strength: 0.6
		
		% Build coupling matrix
		for i=1:N
			x1 = mod(ceil(i/intra_comm_size)-1,n_communities)+1;				% community number
			for j=i:N
				if i~=j									% ignore diagonals
					y1 = mod(ceil(j/intra_comm_size)-1,n_communities)+1;		% community number
					if x1 == y1							% same community
						p = d0/intra_comm_size;
						k = k0;
					else								 % different communities
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
	
end 


%% calculating information atoms

% {
% instantiate variables to store atoms for different coupling matrices and noise correlations
phiid_all_err_coup_mmi = zeros(16, size(coupling_matrices,3), size(error_vec, 2));
phiid_all_err_coup_ccs = zeros(16, size(coupling_matrices,3), size(error_vec, 2));

% instantiate variables to store practical measures for synergistic capacity for different coupling matrices and noise correlations
synergy_capacity_practical = zeros(size(coupling_matrices,3), size(error_vec, 2));
downward_causation_practical = zeros(size(coupling_matrices,3), size(error_vec, 2));
causal_decoupling_practical = zeros(size(coupling_matrices,3), size(error_vec, 2));

% instantiate average covariance matrix
all_average_cov_X = zeros(size(coupling_matrices,3), size(error_vec, 2));
all_average_corr_X = zeros(size(coupling_matrices,3), size(error_vec, 2));

% instantiate matrix for macrio variables

rng(1);
for i = 1:size(coupling_matrices, 3);
	
	if strcmp(network, '8node_mvar_erdoes_renyi') == true
		coupling_matrix = coupling_matrices(:,:,i,:);
		coupling_matrix = squeeze(coupling_matrix(:,:,1,:));
	else coupling_matrix = coupling_matrices(:,:,i);
	end
	% spectral_radius = max(abs(eig(coupling_matrix)));
	disp(i)
	
	for j = 1:length(error_vec)
		
		err = error_vec(j);
		
		if strcmp(network, '8node_mvar_erdoes_renyi') == true;
			phiid_all_err_coup_mmi_temp = zeros(16,50);
			phiid_all_err_coup_ccs_temp = zeros(16,50);
			all_average_cov_X_temp = zeros(1,50);
			all_average_corr_X_temp = zeros(1,50);
			macro_variable = zeros(1, npoints); 
			
			for q = 1:length(fifty_values_per_ER_network);
				X = sim_method(coupling_matrix(:,:,q), npoints, tau, err);
				try
					phiid_all_err_coup_mmi_temp(:,q) = struct2array(PhiIDFull(X, tau, 'MMI'))';
					phiid_all_err_coup_ccs_temp(:,q) = struct2array(PhiIDFull(X, tau, 'ccs'))';
					
					cov_X = cov(X');
					all_average_cov_X_temp(q) = mean(nonzeros(tril(cov_X,-1)), 'all');
					%corr_X = corrcov(cov_X);
					%all_average_corr_X_temp(q) = corrcov(cov_X);
					
					for k = 1:(size(X,1))
						macro_variable = macro_variable + X(k,:);
					end
					
					synergy_capacity_practical(i,j) = EmergencePsi(X', macro_variable');
					downward_causation_practical(i,j) = EmergenceDelta(X', macro_variable');
					causal_decoupling_practical(i,j) = synergy_capacity_practical(i,j) - downward_causation_practical(i,j);
		
				catch 
					phiid_all_err_coup_mmi_temp(:,q) = NaN;
					phiid_all_err_coup_ccs_temp(:,q) = NaN;
					all_average_cov_X_temp(q) = NaN;
				end
			end
			
			% average covariance matrix
			
			all_average_cov_X(i,j) = nanmean(all_average_cov_X_temp);
			all_average_corr_X(i,j) = nanmean(all_average_corr_X_temp);
			phiid_all_err_coup_mmi(:,i,j) = nanmean(real(phiid_all_err_coup_mmi_temp),2);
			phiid_all_err_coup_ccs(:,i,j) = nanmean(real(phiid_all_err_coup_ccs_temp), 2);
		
		elseif strcmp(network, '8node_kuramoto') == true; 	
			[X, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, err, intra_comm_size, n_communities);	
			
			try
				phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
			catch 
				phiid_all_err_coup_mmi(:,i,j) = NaN;
			end
			
			try
				phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
			catch 
				phiid_all_err_coup_ccs(:,i,j) = NaN;
			end
			
			synergy_capacity_practical(i,j) = EmergencePsi(X', sigma_chi');
			downward_causation_practical(i,j) = EmergenceDelta(X', sigma_chi');
			causal_decoupling_practical(i,j) = synergy_capacity_practical(i,j) - downward_causation_practical(i,j);
		
			macro_variables(i,j,:) = sigma_chi;
			synchronies(j,:,:) = synchrony;
			
			% average covariance/correlation matrix
			cov_X = cov(X');
			all_average_cov_X(i,j) = mean(nonzeros(tril(cov_X,-1)), 'all');
			
			corr_X = corrcov(cov_X);
			all_average_corr_X(i,j) = mean(nonzeros(tril(corr_X,-1)), 'all');
			
		elseif strcmp(network, '256node_kuramoto') == true;		% calculating only practical measures with such a big system
			[X, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, err, intra_comm_size, n_communities);
			
			synergy_capacity_practical(i,j) = EmergencePsi(X', sigma_chi');
			downward_causation_practical(i,j) = EmergenceDelta(X', sigma_chi');
			causal_decoupling_practical(i,j) = synergy_capacity_practical(i,j) - downward_causation_practical(i,j);
		
			macro_variables(i,j,:) = sigma_chi;
			synchronies(j,:,:) = synchrony;
			
			% average covariance/correlation matrix
			cov_X = cov(X');
			all_average_cov_X(i,j) = mean(nonzeros(tril(cov_X,-1)), 'all');
			
			corr_X = corrcov(cov_X);
			all_average_corr_X(i,j) = mean(nonzeros(tril(corr_X,-1)), 'all');
		
		else 	
			X = sim_method(coupling_matrix, npoints, tau, err);
		
			% PhiID (synergistic capacity is calculated below)
			phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
			phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		

			% practical measures for causal emergence
		
			% some super simple meaningless macro variable - adding up the micro
			% (for network with community structure, we add the variables of the two communities, respectively, and choose the bigger one for the final macro variable)
		
			macro_variable = zeros(1, npoints);
			if strcmp(network, '8node_mvar_different_architectures') == true;
				if i == 4									% corresponding to network with two community structures
					hub1 = zeros(1, npoints);
					hub2 = zeros(1, npoints);
			
					for k = 1:4
						hub1 = hub1 + X(k,:);
					end  
			
					for k = 5:(size(X,1))
						hub2 = hub2 + X(k,:);
					end
			
					for k = 1:(size(X,2))
						macro_variable(1, k) = max(hub1(1,k), hub2(1,k));
					end
			
				else
					for k = 1:(size(X,1));
						macro_variable = macro_variable + X(k,:);
					end 		
				end 
			else
				for k = 1:(size(X,1));
						macro_variable = macro_variable + X(k,:);
				end 	
			end 
		
			synergy_capacity_practical(i,j) = EmergencePsi(X', macro_variable');
			downward_causation_practical(i,j) = EmergenceDelta(X', macro_variable');
			causal_decoupling_practical(i,j) = synergy_capacity_practical(i,j) - downward_causation_practical(i,j);
		
			% average covariance/correlation matrix
			cov_X = cov(X');
			all_average_cov_X(i,j) = mean(nonzeros(tril(cov_X,-1)), 'all');
			
			corr_X = corrcov(cov_X);
			all_average_corr_X(i,j) = mean(nonzeros(tril(corr_X,-1)), 'all');
	
		end 
		
		
	%spectral_radius
	end 
	
	if ((strcmp(network, '8node_kuramoto') == true) || (strcmp(network, '256node_kuramoto') == true));
		a_string = num2str(coupling_vec(i));
		save([pathout_data network '_synchronies_ '  a_string(3:end) '_' sim_index '.mat'], 'synchronies');
	end
	
end 

% allocating variable names for the atoms in a struct;
% extract single atoms such that rows are the couplings, and columns the errors

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

all_atoms_err_coup_mmi = [];
all_atoms_err_coup_mmi.rtr = squeeze(phiid_all_err_coup_mmi(1,:,:));
all_atoms_err_coup_mmi.rtx = squeeze(phiid_all_err_coup_mmi(2,:,:));
all_atoms_err_coup_mmi.rty = squeeze(phiid_all_err_coup_mmi(3,:,:));
all_atoms_err_coup_mmi.rts = squeeze(phiid_all_err_coup_mmi(4,:,:));
all_atoms_err_coup_mmi.xtr = squeeze(phiid_all_err_coup_mmi(5,:,:));
all_atoms_err_coup_mmi.xtx = squeeze(phiid_all_err_coup_mmi(6,:,:));
all_atoms_err_coup_mmi.xty = squeeze(phiid_all_err_coup_mmi(7,:,:));
all_atoms_err_coup_mmi.xts = squeeze(phiid_all_err_coup_mmi(8,:,:));
all_atoms_err_coup_mmi.ytr = squeeze(phiid_all_err_coup_mmi(9,:,:));
all_atoms_err_coup_mmi.ytx = squeeze(phiid_all_err_coup_mmi(10,:,:));
all_atoms_err_coup_mmi.yty = squeeze(phiid_all_err_coup_mmi(11,:,:));
all_atoms_err_coup_mmi.yts = squeeze(phiid_all_err_coup_mmi(12,:,:));
all_atoms_err_coup_mmi.str = squeeze(phiid_all_err_coup_mmi(13,:,:));
all_atoms_err_coup_mmi.stx = squeeze(phiid_all_err_coup_mmi(14,:,:));
all_atoms_err_coup_mmi.sty = squeeze(phiid_all_err_coup_mmi(15,:,:));
all_atoms_err_coup_mmi.sts = squeeze(phiid_all_err_coup_mmi(16,:,:));

all_atoms_err_coup_ccs = [];
all_atoms_err_coup_ccs.rtr = squeeze(phiid_all_err_coup_ccs(1,:,:));
all_atoms_err_coup_ccs.rtx = squeeze(phiid_all_err_coup_ccs(2,:,:));
all_atoms_err_coup_ccs.rty = squeeze(phiid_all_err_coup_ccs(3,:,:));
all_atoms_err_coup_ccs.rts = squeeze(phiid_all_err_coup_ccs(4,:,:));
all_atoms_err_coup_ccs.xtr = squeeze(phiid_all_err_coup_ccs(5,:,:));
all_atoms_err_coup_ccs.xtx = squeeze(phiid_all_err_coup_ccs(6,:,:));
all_atoms_err_coup_ccs.xty = squeeze(phiid_all_err_coup_ccs(7,:,:));
all_atoms_err_coup_ccs.xts = squeeze(phiid_all_err_coup_ccs(8,:,:));
all_atoms_err_coup_ccs.ytr = squeeze(phiid_all_err_coup_ccs(9,:,:));
all_atoms_err_coup_ccs.ytx = squeeze(phiid_all_err_coup_ccs(10,:,:));
all_atoms_err_coup_ccs.yty = squeeze(phiid_all_err_coup_ccs(11,:,:));
all_atoms_err_coup_ccs.yts = squeeze(phiid_all_err_coup_ccs(12,:,:));
all_atoms_err_coup_ccs.str = squeeze(phiid_all_err_coup_ccs(13,:,:));
all_atoms_err_coup_ccs.stx = squeeze(phiid_all_err_coup_ccs(14,:,:));
all_atoms_err_coup_ccs.sty = squeeze(phiid_all_err_coup_ccs(15,:,:));
all_atoms_err_coup_ccs.sts = squeeze(phiid_all_err_coup_ccs(16,:,:));



save([pathout_data network '_macro_variables_ ' sim_index '.mat'], 'macro_variables');
save([pathout_data network '_all_atoms_err_coup_ccs' sim_index '.mat'], 'all_atoms_err_coup_ccs');
save([pathout_data network '_all_atoms_err_coup_mmi' sim_index '.mat'], 'all_atoms_err_coup_mmi');
save([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
save([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

save([pathout_data network '_all_atoms_err_coup_ccs' sim_index '.mat'], 'all_atoms_err_coup_ccs');
save([pathout_data network '_all_atoms_err_coup_mmi' sim_index '.mat'], 'all_atoms_err_coup_mmi');
save([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
save([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');


%

% allocating variable names for the atoms in a struct;
emergence_practical = [];
emergence_practical = [];

emergence_practical.synergy_capacity_practical = synergy_capacity_practical;
emergence_practical.causal_decoupling_practical = causal_decoupling_practical;
emergence_practical.downward_causation_practical = downward_causation_practical;

save([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');

%}

%% synergistic/emergent capacity, downward causation, causal decoupling

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
synergy_capacity_ccs = all_atoms_err_coup_ccs.str + ...
	all_atoms_err_coup_ccs.stx + all_atoms_err_coup_ccs.sty + all_atoms_err_coup_ccs.sts;

downward_causation_ccs = all_atoms_err_coup_ccs.str + all_atoms_err_coup_ccs.stx + all_atoms_err_coup_ccs.sty;

synergy_capacity_mmi = all_atoms_err_coup_mmi.str + ...
	all_atoms_err_coup_mmi.stx + all_atoms_err_coup_mmi.sty + all_atoms_err_coup_mmi.sts;

downward_causation_mmi = all_atoms_err_coup_mmi.str + all_atoms_err_coup_mmi.stx + all_atoms_err_coup_mmi.sty;

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

save([pathout_plots network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
save([pathout_plots network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
%}

%% plotting 

% axes ticks
if sim_index == '3'
	x_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ... 
		'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
	y_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ... 
		'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
elseif strcmp(network, '2node_mvar') == true;
	x_axis = {'0.09', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
	y_axis = {'0.0045', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '',  '', '0.09',  '', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', ... 
		'0.18', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '',  '', '0.27', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '', ... 
		'0.36', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '0.45', ''};

elseif ((strcmp(network, '8node_mvar_erdoes_renyi') == true) || (strcmp(network, '8node_mvar_global_coupling') == true));
		x_axis = {'0.09', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
		y_axis = {'0.05', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
elseif strcmp(network, '8node_mvar_different_architectures') == true; 
		x_axis = {'0.13', '0.26', '0.39', '0.52', '0.65', '0.78', '0.9'};
		y_axis = {'optimal A', 'optimal B', 'small world', 'two communities', 'fully connected', 'ring', 'uni ring'};

elseif ((strcmp(network, '8node_kuramoto') == true) || (strcmp(network, '256node_kuramoto') == true)); 
		x_axis = {'0.0', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.27', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.54', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.8', ''};
		y_axis = {'0.05', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
end
% double-redundancy & double-synergy

% {
% heatmaps using matlab built-in function

atoms = {all_atoms_err_coup_ccs.rtr, all_atoms_err_coup_mmi.rtr, all_atoms_err_coup_ccs.sts, all_atoms_err_coup_mmi.sts};
file_names = {'_all_err_coup_ccs_rtr', '_all_err_coup_mmi_rtr', '_all_err_coup_ccs_sts', '_all_err_coup_mmi_sts'};
titles = {'double-redundancy ccs', 'double-redundancy mmi', 'double-synergy ccs', 'double-synergy mmi'};

cmin = min([min(atoms{1}(:)), min(atoms{2}(:)), min(atoms{3}(:)), min(atoms{4}(:))]);
cmax = max([max(atoms{1}(:)), max(atoms{2}(:)), max(atoms{3}(:)), max(atoms{4}(:))]);

for i = 1:size(atoms,2)

	figure;

	imagesc(atoms{i});
	colormap(bluewhitered); 
	colorbar;
	
	hColorbar = colorbar;
	set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]))

	xticks(1:100);
	yticks(1:100);
	set(gca,'TickLength',[0 0])
	yticklabels(y_axis);
	xticklabels(x_axis);
	
	if sim_index == '3';
		ylabel = 'zero coupling';
		xlabel = 'zero noise correlation';
	
	else 
		xlabel = 'noise correlation';
		if strcmp(network, '2node') == true;
			ylabel = 'coupling strength';
		elseif strcmp(network, '8node_erdoes_renyi') == true;
			ylabel = 'density';
		elseif strcmp(network, '8node_global_coupling') == true;
			ylabel = 'global coupling factor';
		else 
			ylabel = 'network architecture';
		end
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end

	

%}

% synergistic capacity, downward causation, causal decoupling

% {
% heatmaps using matlab built-in function:

if strcmp(network, '256node_kuramoto') == true;
	atoms = {synergy_capacity_practical, downward_causation_practical, causal_decoupling_practical, all_average_cov_X, all_average_corr_X};
	file_names = {'_all_err_coup_practical_synergy_capacity', '_all_err_coup_practical_downward_causation', 
		'_all_err_coup_practical_causal_decoupling', '_all_err_coup_average_cov_X', '_all_err_coup_average_corr_X'};
	titles = {'synergy capacity ccs', 'synergy capacity mmi', 'synergy capacity practical', 'downward causation ccs', 
		'downward causation mmi', 'downward causation practical', 'causal decoupling ccs', 'causal decoupling mmi', 
		'causal decoupling practical', 'average covariance X', 'average correlation X'};
else
	atoms = {synergy_capacity_ccs, synergy_capacity_mmi, synergy_capacity_practical, downward_causation_ccs, downward_causation_mmi, downward_causation_practical, causal_decoupling_ccs, causal_decoupling_mmi, causal_decoupling_practical, all_average_cov_X, all_average_corr_X};
	file_names = {'_all_err_coup_ccs_synergy_capacity', '_all_err_coup_mmi_synergy_capacity',  '_all_err_coup_practical_synergy_capacity', 
		'_all_err_coup_ccs_downward_causation', '_all_err_coup_mmi_downward_causation', '_all_err_coup_practical_downward_causation', 
		'_all_err_coup_ccs_causal_decoupling', '_all_err_coup_mmi_causal_decoupling', '_all_err_coup_practical_causal_decoupling', 
		'_all_err_coup_average_cov_X', '_all_err_coup_average_corr_X'};
	titles = {'synergy capacity ccs', 'synergy capacity mmi', 'synergy capacity practical', 'downward causation ccs', 'downward causation mmi', 
		'downward causation practical', 'causal decoupling ccs', 'causal decoupling mmi', 'causal decoupling practical', 'average covariance X', 
		'average correlation X'};
end 


for i = 1:size(atoms,2)
	
	figure;

	imagesc(atoms{i});
	colormap(bluewhitered); 
	colorbar;
	
	hColorbar = colorbar;
	set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]))

	xticks(1:100);
	yticks(1:100);
	set(gca,'TickLength',[0 0])
	yticklabels(y_axis);
	xticklabels(x_axis);
	
	if sim_index == '3'
		ylabel = 'zero coupling';
		xlabel = 'zero noise correlation';
	
	else 
		xlabel = 'noise correlation';
		if strcmp(network, '2node_mvar') == true; 
			ylabel = 'coupling strength';
		elseif strcmp(network, '8node_mvar_erdoes_renyi') == true;
			ylabel = 'density';
		elseif strcmp(network, '8node_mvar_global_coupling') == true;
			ylabel = 'global coupling factor';
		else 
			ylabel = 'network architecture';
		end
	end 
	
	title(titles{i});
	exportgraphics(gcf, [PATHOUT2 network file_names{i} sim_index '.png']);

end
%}

close all;


%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 8-node kuramoto oscillators, with fixed A, and varying beta

if strcmp(network, '8node_kuramoto') == true;
	
clear ylabel;
clear xlabel;

load([pathout_data network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
load([pathout_data network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');

% coupling_vec(20) = 0.2131
A = [5, 20, 40, 60, 80, 95];

for p = 1:length(A);

g = A(p);

blubb1 = emergence_ccs.synergy_capacity_ccs(g,:);
figure;
scatter(error_vec, blubb1, 60, 'filled');
title(['emergence capcity ccs, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity ccs');
xlabel('beta');

a_string = num2str(coupling_vec(g));
a_string = a_string(3:end);
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_ccs_' a_string '_' sim_index '.png']);

blubb2 = emergence_mmi.synergy_capacity_mmi(g,:);
figure;
scatter(error_vec, blubb2, 60, 'filled');
title(['emergence capacity mmi, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity mmi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_mmi_' a_string '_' sim_index '.png']);

blubb3 = emergence_practical.synergy_capacity_practical(g,:);
figure;
scatter(error_vec, blubb3, 60, 'filled');
title(['emergence capacity practical, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity practical');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_practical_' a_string '_' sim_index '.png']);

load([pathout_data network '_synchronies_ '  a_string '_' sim_index '.mat'], 'synchronies');
%load([PATHOUT3 network '_synchronies_ '  '21313' '_' sim_index '.mat'], 'synchronies');

% metastability
sigma_met_mean = [];
for i = 1:length(error_vec);
	sigma_met = squeeze(synchronies(i,:,:));
	sigma_met_mean(i) = mean(var(sigma_met'));
end

% chimera states
sigma_chi_mean = [];
for i = 1:length(error_vec);
	sigma_chi = squeeze(synchronies(i,:,:));
	sigma_chi_mean(i) = mean(var(sigma_chi));
end

figure;
scatter(error_vec, sigma_chi_mean, 60, 'filled');
title(['sigma chi mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma chi mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_sigma_chi_mean_' a_string '_' sim_index '.png']);

figure;
scatter(error_vec, sigma_met_mean, 60, 'filled');
title(['sigma met mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma met mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_sigma_met_mean_' a_string '_' sim_index '.png']);

close all;
end 

end 
%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 256-node kuramoto oscillators, with fixed A, and varying beta

if strcmp(network, '256node_kuramoto') == true;
	
clear ylabel;
clear xlabel;

load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');

% coupling_vec(20) = 0.2131
A = [5, 20, 40, 60, 80, 95];

for p = 1:length(A);

g = A(p);

blubb3 = emergence_practical.synergy_capacity_practical(g,:);
figure;
scatter(error_vec, blubb3, 60, 'filled');
title(['emergence capacity practical, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity practical');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_practical_' a_string '_' sim_index '.png']);

load([pathout_data network '_synchronies_ '  a_string '_' sim_index '.mat'], 'synchronies');
%load([PATHOUT3 network '_synchronies_ '  '21313' '_' sim_index '.mat'], 'synchronies');

% metastability
sigma_met_mean = [];
for i = 1:length(error_vec);
	sigma_met = squeeze(synchronies(i,:,:));
	sigma_met_mean(i) = mean(var(sigma_met'));
end

% chimera states
sigma_chi_mean = [];
for i = 1:length(error_vec);
	sigma_chi = squeeze(synchronies(i,:,:));
	sigma_chi_mean(i) = mean(var(sigma_chi));
end

figure;
scatter(error_vec, sigma_chi_mean, 60, 'filled');
title(['sigma chi mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma chi mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_sigma_chi_mean_' a_string '_' sim_index '.png']);

figure;
scatter(error_vec, sigma_met_mean, 60, 'filled');
title(['sigma met mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma met mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_sigma_met_mean_' a_string '_' sim_index '.png']);

close all;
end 

end


%% generating built-in heatmaps with colormap parula

%{
for i = 1:size(atoms,2)
	
	figure;
	clf
	h = heatmap(atoms{i}, 'Colormap', parula, 'ColorbarVisible', 'on', 'CellLabelColor', 'none') ;
	h.YDisplayLabels = y_axis;
	h.XDisplayLabels = x_axis;
	%caxis([0, 0.1]);
	
	if sim_index == '3'
		h.YLabel = 'zero coupling';
		h.XLabel = 'zero noise correlation';
	
	else 
		h.XLabel = 'noise correlation';
		if strcmp(network, '2node') == true; 
			h.YLabel = 'coupling strength';
		elseif strcmp(network, '8node_erdoes_renyi') == true;
			h.YLabel = 'density';
		elseif strcmp(network, '8node_erdoes_global_coupling') == true;
			h.YLabel = 'global coupling factor';
		else 
			h.YLabel = 'network architecture';
		end
	end 
	
	title(titles{i});
	exportgraphics(gcf, [PATHOUT2 network file_names{i} sim_index '.png']);

end
%}