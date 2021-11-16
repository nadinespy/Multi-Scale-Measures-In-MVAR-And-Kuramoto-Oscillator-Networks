%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?
% - add integrated information measures?

%% 8-NODE MULTIVARIATE AUTOREGRESSIVE (MVAR) ERDOES-RENYI TIME-SERIES WITH DIFFERENT DENSITIES & MACRO VARIABLES

% This script implements synergy capacity for 8-node MVAR Erdoes-Renyi models with different densities and noise correlations, and two different macro variables 
% (one time being a summation of the raw values of X, and one time being a summation of the exponential of X).

clear all;
clear java;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

pathout_data = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/8node_mvar_erdoes_renyi/';
pathout_plots = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/8node_mvar_erdoes_renyi/';

%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
npoints = 2000;
tau = 1;

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_random(), chimera_metastable_model())
sim_method = @statdata_coup_errors1;

% save plots & matrices according to simulation method (options: '1' (for statdata_coup_errors1()), '2' (for statdata_coup_errors2()), '3' (for statdata_random(), '4' (for metastable_chimera_model())
sim_index = '1';

% network (options: '2node_mvar' for 2-node network with 100 different coupling strengths & noise correlations (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations (if choosing sim_index = 3);
% '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations (if choosing sim_index = 3));
% '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations (if choosing sim_index = 1 or 2)
% '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations (if choosing sim_index = 1 or 2)
% '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 4) 
% '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 5) 

network = '8node_mvar_erdoes_renyi';

%% load files (if already existent, to, e. g., only create plots)

%{

load([pathout_data network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
load([pathout_data network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');
load([pathout_data network '_all_atoms_err_coup_ccs' sim_index '.mat'], 'all_atoms_err_coup_ccs');
load([pathout_data network '_all_atoms_err_coup_mmi' sim_index '.mat'], 'all_atoms_err_coup_mmi');
load([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
load([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

synergy_capacity_mmi = emergence_mmi.synergy_capacity_mmi;
downward_causation_mmi = emergence_mmi.downward_causation_mmi;
causal_decoupling_mmi = emergence_mmi.causal_decoupling_mmi;

synergy_capacity_ccs = emergence_ccs.synergy_capacity_ccs;
downward_causation_ccs = emergence_ccs.downward_causation_ccs;
causal_decoupling_ccs = emergence_ccs.causal_decoupling_ccs;

synergy_capacity_practical_linear = emergence_practical.synergy_capacity_practical_linear; 
downward_causation_practical_linear = emergence_practical.downward_causation_practical_linear; 
causal_decoupling_practical_linear = emergence_practical.causal_decoupling_practical_linear;

synergy_capacity_practical_exponential = emergence_practical.synergy_capacity_practical_exponential; 
downward_causation_practical_exponential = emergence_practical.downward_causation_practical_exponential; 
causal_decoupling_practical_exponential = emergence_practical.causal_decoupling_practical_exponential;

%}

%% create density matrices & noise correlation vectors

% {

density_vec = linspace(0.05,0.9, 10);
error_vec = linspace(0.01, 0.9, 10);
nb_ER_networks = 50;
	
density_matrices = [];
	
for n = 1:length(error_vec);						
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
			density_matrices(:,:,o,p) = ER_network_normalized;
			
		end
	end
end


%% calculating information atoms & practical measures

% {

% instantiate variables to store practical measures for synergistic capacity for different coupling matrices and noise correlations

% macro variable: summation of exponential values of X
synergy_capacity_practical_exponential = zeros(size(density_matrices,3), size(error_vec, 2));
downward_causation_practical_exponential = zeros(size(density_matrices,3), size(error_vec, 2));
causal_decoupling_practical_exponential = zeros(size(density_matrices,3), size(error_vec, 2));

% macro variable: summation of raw values of X
synergy_capacity_practical_linear = zeros(size(density_matrices,3), size(error_vec, 2));
downward_causation_practical_linear = zeros(size(density_matrices,3), size(error_vec, 2));
causal_decoupling_practical_linear = zeros(size(density_matrices,3), size(error_vec, 2));

% instantiate average covariance/correlation matrix
all_average_cov_X = zeros(size(density_matrices,3), size(error_vec, 2));
all_average_corr_X = zeros(size(density_matrices,3), size(error_vec, 2));

rng(1);
for i = 1:size(density_matrices, 3);
	
	density_matrix = density_matrices(:,:,i,:);
	density_matrix = squeeze(density_matrix(:,:,1,:));
	disp(i)
	
	for j = 1:length(error_vec)
		
		err = error_vec(j);

		phiid_all_err_coup_mmi_temp = zeros(16,50);
		phiid_all_err_coup_ccs_temp = zeros(16,50);
		all_average_cov_X_temp = zeros(1,50);
		all_average_corr_X_temp = zeros(1,50);
		synergy_capacity_practical_linear_temp = zeros(1,50);
		downward_causation_practical_linear_temp = zeros(1,50);
		causal_decoupling_practical_linear_temp = zeros(1,50);
		synergy_capacity_practical_exponential_temp = zeros(1,50);
		downward_causation_practical_exponential_temp = zeros(1,50);
		causal_decoupling_practical_exponential_temp = zeros(1,50);
			
		for q = 1:length(fifty_values_per_ER_network);
			X = sim_method(density_matrix(:,:,q), npoints, tau, err);
			
			macro_variable_linear_temp = zeros(1, npoints); 
			macro_variable_exponential_temp = zeros(1, npoints); 
			
			try
				phiid_all_err_coup_mmi_temp(:,q) = struct2array(PhiIDFull(X, tau, 'MMI'))';
				phiid_all_err_coup_ccs_temp(:,q) = struct2array(PhiIDFull(X, tau, 'ccs'))';
					
				cov_X = cov(X');
				all_average_cov_X_temp(q) = mean(nonzeros(tril(cov_X,-1)), 'all');
				corr_X = corrcov(cov_X);
				all_average_corr_X_temp(q) = mean(nonzeros(tril(corr_X,-1)), 'all');
					
				for k = 1:(size(X,1))
					macro_variable_linear_temp = macro_variable_linear_temp + X(k,:);
				end
					
				synergy_capacity_practical_linear_temp(q) = EmergencePsi(X', macro_variable_linear_temp');
				downward_causation_practical_linear_temp(q) = EmergenceDelta(X', macro_variable_linear_temp');
				causal_decoupling_practical_linear_temp(q) = synergy_capacity_practical_linear_temp(q) - downward_causation_practical_linear_temp(q);
				
				for k = 1:(size(X,1))
					macro_variable_exponential_temp = macro_variable_exponential_temp + exp(X(k,:));
				end
					
				synergy_capacity_practical_exponential_temp(q) = EmergencePsi(X', macro_variable_exponential_temp');
				downward_causation_practical_exponential_temp(q) = EmergenceDelta(X', macro_variable_exponential_temp');
				causal_decoupling_practical_exponential_temp(q) = synergy_capacity_practical_exponential_temp(q) - downward_causation_practical_exponential_temp(q);
		
			catch 
				phiid_all_err_coup_mmi_temp(:,q) = NaN;
				phiid_all_err_coup_ccs_temp(:,q) = NaN;
				all_average_cov_X_temp(q) = NaN;
				all_average_corr_X_temp(q) = NaN;
				synergy_capacity_practical_exponential_temp(q) = NaN;
				downward_causation_practical_exponential_temp(q) = NaN;
				causal_decoupling_practical_exponential_temp(q) = NaN;
				synergy_capacity_practical_linear_temp(q) = NaN;
				downward_causation_practical_linear_temp(q) = NaN;
				causal_decoupling_practical_linear_temp(q) = NaN;

			end
		end 
			
		% averages of all 50 ER networks
		all_average_cov_X(i,j) = nanmean(all_average_cov_X_temp);
		all_average_corr_X(i,j) = nanmean(all_average_corr_X_temp);
		phiid_all_err_coup_mmi(:,i,j) = nanmean(real(phiid_all_err_coup_mmi_temp),2);
		phiid_all_err_coup_ccs(:,i,j) = nanmean(real(phiid_all_err_coup_ccs_temp), 2);
		synergy_capacity_practical_exponential(i,j) = nanmean(synergy_capacity_practical_exponential_temp);
		downward_causation_practical_exponential(i,j) = nanmean(downward_causation_practical_exponential_temp);
		causal_decoupling_practical_exponential(i,j) = nanmean(causal_decoupling_practical_exponential_temp);
		synergy_capacity_practical_linear(i,j) = nanmean(synergy_capacity_practical_linear_temp);
		downward_causation_practical_linear(i,j) = nanmean(downward_causation_practical_linear_temp);
		causal_decoupling_practical_linear(i,j) = nanmean(causal_decoupling_practical_linear_temp);
		
	end
end 

save([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
save([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');

%% storing information atoms & practical measures for different macro variables in struct files

% practical measures for different macro variables
emergence_practical = [];

emergence_practical.synergy_capacity_practical_linear = synergy_capacity_practical_linear;
emergence_practical.causal_decoupling_practical_linear = causal_decoupling_practical_linear;
emergence_practical.downward_causation_practical_linear = downward_causation_practical_linear;
emergence_practical.synergy_capacity_practical_exponential = synergy_capacity_practical_exponential;
emergence_practical.causal_decoupling_practical_exponential = causal_decoupling_practical_exponential;
emergence_practical.downward_causation_practical_exponential = downward_causation_practical_exponential;

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

save([pathout_data network '_all_atoms_err_coup_ccs' sim_index '.mat'], 'all_atoms_err_coup_ccs');
save([pathout_data network '_all_atoms_err_coup_mmi' sim_index '.mat'], 'all_atoms_err_coup_mmi');
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
	x_axis = {'0.09', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
	y_axis = {'0.05', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
		'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
end

% double-redundancy & double-synergy

% {
% heatmaps using imagesc()

atoms = {all_atoms_err_coup_ccs.rtr, ...
	all_atoms_err_coup_mmi.rtr, ...
	all_atoms_err_coup_ccs.sts, ...
	all_atoms_err_coup_mmi.sts};
file_names = {'_all_err_coup_ccs_rtr', ...
	'_all_err_coup_mmi_rtr', ...
	'_all_err_coup_ccs_sts', ...
	'_all_err_coup_mmi_sts'};
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
		xlabel('noise correlation');
		ylabel('density');
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end

%}

% synergistic capacity, downward causation, causal decoupling

% {
% heatmaps using imagesc()

atoms = {synergy_capacity_ccs, synergy_capacity_mmi, ...
	synergy_capacity_practical_linear, ...
	synergy_capacity_practical_exponential, ...
	downward_causation_ccs, downward_causation_mmi, ...
	downward_causation_practical_linear, ...
	downward_causation_practical_exponential, ...
	causal_decoupling_ccs, causal_decoupling_mmi, ...
	causal_decoupling_practical_linear, ...
	causal_decoupling_practical_exponential, ...
	all_average_cov_X, ...
	all_average_corr_X};

file_names = {'_all_err_coup_ccs_synergy_capacity', ...
	'_all_err_coup_mmi_synergy_capacity',  ...
	'_all_err_coup_synergy_capacity_practical_linear', ...
	'_all_err_coup_synergy_capacity_practical_exponential', ...
	'_all_err_coup_ccs_downward_causation', ...
	'_all_err_coup_mmi_downward_causation', ...
	'_all_err_coup_downward_causation_practical_linear', ...
	'_all_err_coup_downward_causation_practical_exponential', ...
	'_all_err_coup_ccs_causal_decoupling', ...
	'_all_err_coup_mmi_causal_decoupling', ...
	'_all_err_coup_causal_decoupling_practical_linear', ...
	'_all_err_coup_causal_decoupling_practical_exponential', ...
	'_all_err_coup_average_cov_X', ...
	'_all_err_coup_average_corr_X'};

titles = {'synergy capacity ccs', 'synergy capacity mmi', ...
	'synergy capacity practical linear', ...
	'synergy capacity practical exponential', ...
	'downward causation ccs', 'downward causation mmi', ...
	'downward causation practical linear', ...
	'downward causation practical exponential', ...
	'causal decoupling ccs', ...
	'causal decoupling mmi', ...
	'causal decoupling practical linear', ...
	'causal decoupling practical exponential', ...
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
		xlabel('noise correlation');
		ylabel('density');
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end
%}

%close all;
