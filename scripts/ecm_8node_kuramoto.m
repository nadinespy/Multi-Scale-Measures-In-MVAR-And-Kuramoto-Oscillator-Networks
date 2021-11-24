%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - separate calculations for phiid-based and practical measures 
% - fill struct file in a loop?
% - add integrated information measures?

%% 8-NODE KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS & MACRO VARIABLES

% This script implements synergy capacity for 8-node Kuramoto oscillators with different couplings and phase lags, and two different macro variables 
% (variance of synchronies & global average pairwise synchrony between communities).

clear all;
clear java;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

pathout_data = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_kuramoto/';
pathout_plots = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_kuramoto/';


%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
npoints = 2000;
tau = 1;

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_random(), chimera_metastable_model())
sim_method = @chimera_metastable_model;

% save plots & matrices according to simulation method (options: '1' (for statdata_coup_errors1()), '2' (for statdata_coup_errors2()), '3' (for statdata_random(), '4' (for metastable_chimera_model())
sim_index = '4';

% network (options: '2node_mvar' for 2-node network with 100 different coupling strengths & noise correlations (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations (if choosing sim_index = 3);
% '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations (if choosing sim_index = 3));
% '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations (if choosing sim_index = 1 or 2)
% '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations (if choosing sim_index = 1 or 2)
% '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 4) 
% '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 5) 

network = '8node_kuramoto';

%% load files (if already existent, to, e. g., only create plots)

% {

load([pathout_data network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
load([pathout_data network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');
load([pathout_data network '_all_atoms_beta_coup_ccs' sim_index '.mat'], 'all_atoms_beta_coup_ccs');
load([pathout_data network '_all_atoms_beta_coup_mmi' sim_index '.mat'], 'all_atoms_beta_coup_mmi');
load([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
load([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

synergy_capacity_mmi = emergence_mmi.synergy_capacity_mmi;
downward_causation_mmi = emergence_mmi.downward_causation_mmi;
causal_decoupling_mmi = emergence_mmi.causal_decoupling_mmi;

synergy_capacity_ccs = emergence_ccs.synergy_capacity_ccs;
downward_causation_ccs = emergence_ccs.downward_causation_ccs;
causal_decoupling_ccs = emergence_ccs.causal_decoupling_ccs;

synergy_capacity_practical_sigma_chi = emergence_practical.synergy_capacity_practical_sigma_chi; 
downward_causation_practical_sigma_chi = emergence_practical.downward_causation_practical_sigma_chi; 
causal_decoupling_practical_sigma_chi = emergence_practical.causal_decoupling_practical_sigma_chi;

synergy_capacity_practical_pairwise_synchrony = emergence_practical.synergy_capacity_practical_pairwise_synchrony; 
downward_causation_practical_pairwise_synchrony = emergence_practical.downward_causation_practical_pairwise_synchrony; 
causal_decoupling_practical_pairwise_synchrony = emergence_practical.causal_decoupling_practical_pairwise_synchrony;

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
	
synchronies = zeros(length(beta_vec), n_communities, npoints);
	
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

% instantiate variables to store practical measures for synergistic capacity for different coupling matrices and noise correlations

% macro variable: variance of synchronies 
synergy_capacity_practical_sigma_chi = zeros(size(coupling_matrices,3), size(beta_vec, 2));
downward_causation_practical_sigma_chi = zeros(size(coupling_matrices,3), size(beta_vec, 2));
causal_decoupling_practical_sigma_chi = zeros(size(coupling_matrices,3), size(beta_vec, 2));

% macro variable: global average pairwise synchrony between communities
synergy_capacity_practical_pairwise_synchrony = zeros(size(coupling_matrices,3), size(beta_vec, 2));
downward_causation_practical_pairwise_synchrony = zeros(size(coupling_matrices,3), size(beta_vec, 2));
causal_decoupling_practical_pairwise_synchrony = zeros(size(coupling_matrices,3), size(beta_vec, 2));

% average covariance/correlation matrix
all_average_cov_X = zeros(size(coupling_matrices,3), size(beta_vec, 2));
all_average_corr_X = zeros(size(coupling_matrices,3), size(beta_vec, 2));

rng(1);
for i = 1:size(coupling_matrices, 3);
	
	coupling_matrix = coupling_matrices(:,:,i);
	disp(i)
	
	for j = 1:length(beta_vec)
		
		beta = beta_vec(j);
		
		[X, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);	
		
		% PhiID
		try
			phiid_all_beta_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
		catch 
			phiid_all_beta_coup_mmi(:,i,j) = NaN;
		end
			
		try
			phiid_all_beta_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		catch 
			phiid_all_beta_coup_ccs(:,i,j) = NaN;
		end
			
		synchronies(j,:,:) = synchrony;				% rows: betas, columns: communities, 3rd dimension: time-points
		
		% practical measures for causal emergence: variance of synchronies & global average pairwise synchrony between communities
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
		
		synergy_capacity_practical_sigma_chi(i,j) = EmergencePsi(X', sigma_chi');
		synergy_capacity_practical_pairwise_synchrony(i,j) = EmergencePsi(X', grand_mean_pairwise_synchrony');
		downward_causation_practical_sigma_chi(i,j) = EmergenceDelta(X', sigma_chi');
		downward_causation_practical_pairwise_synchrony(i,j) = EmergenceDelta(X', grand_mean_pairwise_synchrony');
		causal_decoupling_practical_sigma_chi(i,j) = synergy_capacity_practical_sigma_chi(i,j) - downward_causation_practical_sigma_chi(i,j);
		causal_decoupling_practical_pairwise_synchrony(i,j) = synergy_capacity_practical_pairwise_synchrony(i,j) - downward_causation_practical_pairwise_synchrony(i,j);

		% average covariance/correlation matrix
		cov_X = cov(X');
		all_average_cov_X(i,j) = mean(nonzeros(tril(cov_X,-1)), 'all');
			
		corr_X = corrcov(cov_X);
		all_average_corr_X(i,j) = mean(nonzeros(tril(corr_X,-1)), 'all');
			
	end
	
	a_string = num2str(coupling_vec(i));
	save([pathout_data network '_synchronies_ '  a_string(3:end) '_' sim_index '.mat'], 'synchronies');
	
end 

save([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
save([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

%% storing information atoms & practical measures for different macro variables in struct files

% practical measures for different macro variables

emergence_practical = [];

emergence_practical.synergy_capacity_practical_sigma_chi = synergy_capacity_practical_sigma_chi;
emergence_practical.causal_decoupling_practical_sigma_chi = causal_decoupling_practical_sigma_chi;
emergence_practical.downward_causation_practical_sigma_chi = downward_causation_practical_sigma_chi;
emergence_practical.synergy_capacity_practical_pairwise_synchrony = synergy_capacity_practical_pairwise_synchrony;
emergence_practical.causal_decoupling_practical_pairwise_synchrony = causal_decoupling_practical_pairwise_synchrony;
emergence_practical.downward_causation_practical_pairwise_synchrony = downward_causation_practical_pairwise_synchrony;

save([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');

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

save([pathout_data network '_all_atoms_beta_coup_ccs' sim_index '.mat'], 'all_atoms_beta_coup_ccs');
save([pathout_data network '_all_atoms_beta_coup_mmi' sim_index '.mat'], 'all_atoms_beta_coup_mmi');

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

save([pathout_data network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
save([pathout_data network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
%}

%% plotting

% axes ticks
if sim_index == '3'
	x_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ... 
		'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
	y_axis = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''  ... 
		'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' };
else
	x_axis = {'0.0', '', '0.18', '', '0.36', '', '0.55', '', '', '0.8'};
	y_axis = {'0.0', '', '0.24', '', '0.43', '', '0.62', '', '', '0.9'};
end

% double-redundancy & double-synergy

% {
% heatmaps using imagesc()

atoms = {all_atoms_beta_coup_ccs.rtr, ...
	all_atoms_beta_coup_mmi.rtr, ...
	all_atoms_beta_coup_ccs.sts, ...
	all_atoms_beta_coup_mmi.sts};
file_names = {'_all_beta_coup_ccs_rtr', ...
	'_all_beta_coup_mmi_rtr', ...
	'_all_beta_coup_ccs_sts', ...
	'_all_beta_coup_mmi_sts'};
titles = {'double-redundancy ccs', ...
	'double-redundancy mmi', ...
	'double-synergy ccs', ...
	'double-synergy mmi'};

% cmin = min([min(atoms{1}(:)), min(atoms{2}(:)), min(atoms{3}(:)), min(atoms{4}(:))]);
% cmax = max([max(atoms{1}(:)), max(atoms{2}(:)), max(atoms{3}(:)), max(atoms{4}(:))]);

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
	
	if sim_index == '3';
		ylabel('zero coupling');
		xlabel('zero noise correlation');
	else 
		ylabel('A');
		xlabel('beta');
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end


%}

% synergistic capacity, downward causation, causal decoupling

% {
% heatmaps using matlab built-in function:

atoms = {synergy_capacity_ccs, ...
	synergy_capacity_mmi, ...
	synergy_capacity_practical_pairwise_synchrony, ...
	synergy_capacity_practical_sigma_chi, ...
	downward_causation_ccs, ...
	downward_causation_mmi, ...
	downward_causation_practical_pairwise_synchrony, ...
	downward_causation_practical_sigma_chi, ...
	causal_decoupling_ccs, ...
	causal_decoupling_mmi, ...
	causal_decoupling_practical_pairwise_synchrony, ...
	causal_decoupling_practical_sigma_chi, ...
	all_average_cov_X, ...
	all_average_corr_X};

file_names = {'_all_beta_coup_ccs_synergy_capacity', ...
	'_all_beta_coup_mmi_synergy_capacity', ...
	'_all_beta_coup_practical_synergy_capacity_pairwise_synchrony', ...
	'_all_beta_coup_practical_synergy_capacity_sigma_chi', ...
	'_all_beta_coup_ccs_downward_causation', ...
	'_all_beta_coup_mmi_downward_causation', ...
	'_all_beta_coup_practical_downward_causation_pairwise_synchrony', ...
	'_all_beta_coup_practical_downward_causation_sigma_chi', ...
	'_all_beta_coup_ccs_causal_decoupling', ...
	'_all_beta_coup_mmi_causal_decoupling', ...
	'_all_beta_coup_practical_causal_decoupling_pairwise_synchrony', ...
	'_all_beta_coup_practical_causal_decoupling_sigma_chi', ...
	'_all_beta_coup_average_cov_X', ...
	'_all_beta_coup_average_corr_X'};

titles = {'synergy capacity ccs', ...
	'synergy capacity mmi', ...
	'synergy capacity practical pairwise synchrony', ...
	'synergy capacity practical sigma chi', ...
	'downward causation ccs', ...
	'downward causation mmi', ...
	'downward causation practical pairwise synchrony', ...
	'downward causation practical sigma chi', ...
	'causal decoupling ccs', ...
	'causal decoupling mmi', ...
	'causal decoupling practical pairwise synchrony', ...
	'causal decoupling practical sigma chi', ...
	'average covariance X', ...
	'average correlation X'};

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
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end
%}

%close all;


%% scatter plots for emergence capacity, sigma met mean & sigma chi mean in 8-node kuramoto oscillators, with fixed A, and varying beta

A = [2, 4, 6, 8, 10];

for p = 1:length(A);

g = A(p);

synergy_capacity_ccs_temp = emergence_ccs.synergy_capacity_ccs(g,:);
figure;
scatter(beta_vec, synergy_capacity_ccs_temp , 60, 'filled');
title(['emergence capcity ccs, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity ccs');
xlabel('beta');

a_string = num2str(coupling_vec(g));
a_string = a_string(3:end);
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_ccs_' a_string '_' sim_index '.png']);

synergy_capacity_mmi_temp = emergence_mmi.synergy_capacity_mmi(g,:);
figure;
scatter(beta_vec, synergy_capacity_mmi_temp, 60, 'filled');
title(['emergence capacity mmi, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity mmi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_mmi_' a_string '_' sim_index '.png']);

synergy_capacity_practical_pairwise_synchrony_temp = emergence_practical.synergy_capacity_practical_pairwise_synchrony(g,:);
figure;
scatter(beta_vec, synergy_capacity_practical_pairwise_synchrony_temp, 60, 'filled');
title(['emergence capacity practical, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity practical pairwise synchrony');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_practical_pairwise_synchrony_' a_string '_' sim_index '.png']);

synergy_capacity_practical_sigma_chi_temp = emergence_practical.synergy_capacity_practical_sigma_chi(g,:);
figure;
scatter(beta_vec, synergy_capacity_practical_sigma_chi_temp, 60, 'filled');
title(['emergence capacity practical sigma chi, A = ' num2str(coupling_vec(g))]);
ylabel('emergence capacity practical sigma chi');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_synergy_capacity_practical_sigma_chi_' a_string '_' sim_index '.png']);

load([pathout_data network '_synchronies_ '  a_string '_' sim_index '.mat'], 'synchronies');
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
exportgraphics(gcf, [pathout_plots network '_sigma_chi_mean_' a_string '_' sim_index '.png']);

figure;
scatter(beta_vec, sigma_met_mean, 60, 'filled');
title(['sigma met mean, A = ' num2str(coupling_vec(g))]);
ylabel('sigma met mean');
xlabel('beta');
exportgraphics(gcf, [pathout_plots network '_sigma_met_mean_' a_string '_' sim_index '.png']);

end 

close all;
