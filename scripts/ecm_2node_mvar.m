%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling matrix in the saved mat-file?
% - fill struct file in a loop?
% - add integrated information measures?

%% 2-NODE MULTIVARIATE AUTOREGRESSIVE TIME-SERIES (MVAR)

% This script implements synergy capacity for 2-node MVAR models with differing connection strengths and noise correlations, and two different macro variables 
% (one time being a summation of the raw values of X, and one time being a summation of the exponential of X).

clear all;
clear java;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/ReconcilingEmergences-master'
javaaddpath('infodynamics.jar');

pathout_data = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/analyses/2node_mvar/';
pathout_plots = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/results/plots/2node_mvar/';

%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
npoints = 2000;
tau = 3;

% simulation method (options: statdata_coup_errors1(), statdata_coup_errors2(), statdata_coup_errors3(), statdata_random(), chimera_metastable_model())
sim_method = @statdata_coup_errors1;

% save plots & matrices according to simulation method (options: '1' (for statdata_coup_errors1()), '2' (for statdata_coup_errors2()), '3' (for statdata_random(), '4' (for metastable_chimera_model(), '5' for statdata_coup_errors3())
sim_index = '1';

% network (options: '2node_mvar' for 2-node network with 100 different coupling strengths & noise correlations (if choosing sim_index = 1 or 2) OR random 2-node network with 100 zero couplings & 100 zero noise correlations (if choosing sim_index = 3);
% '8node_mvar_different_architectures' for 8-node networks with 6 different architectures & noise correlations (if choosing sim_index = 1 or 2) OR random 8-node networks with 100 zero couplings & 100 zero correlations (if choosing sim_index = 3));
% '8node_mvar_erdoes_renyi' for 8-node Erdös-Renyi networks with 100 different densities & noise correlations (if choosing sim_index = 1 or 2)
% '8node_mvar_global_coupling' for phi-optimal network with 100 different global coupling factors & noise correlations (if choosing sim_index = 1 or 2)
% '8node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 4) 
% '256node_kuramoto' for metastable chimera states with 100 different intra- and intercommunity coupling strengths & betas (if choosing sim_index = 5) 

network = '2node_mvar';
method_DD = 'Gaussian';
tau_steps = 1;

%% load files (if already existent, to, e. g., only create plots)

%{
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

	synergy_capacity_practical_linear = emergence_practical.synergy_capacity_practical_linear; 
	downward_causation_practical_linear = emergence_practical.downward_causation_practical_linear; 
	causal_decoupling_practical_linear = emergence_practical.causal_decoupling_practical_linear;

	synergy_capacity_practical_exponential = emergence_practical.synergy_capacity_practical_exponential; 
	downward_causation_practical_exponential = emergence_practical.downward_causation_practical_exponential; 
	causal_decoupling_practical_exponential = emergence_practical.causal_decoupling_practical_exponential;

%}

%% create coupling matrices & noise correlation vectors

% {
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

%}

%% calculating information atoms & practical measures

% {

% instantiate variables to store practical measures for synergistic capacity for different coupling matrices and noise correlations

nbins = 100;

% calculate PhiID & practical measures
rng(1);
for i = 1:size(coupling_matrices, 3);
	
	coupling_matrix = coupling_matrices(:,:,i);
	disp(i)
	
	for j = 1:length(error_vec)

		err = error_vec(j);

		X = sim_method(coupling_matrix, npoints, tau, err);

		% PhiID 
		% {
		phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
		phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		%}
		
		% practical measures for causal emergence - some super simple meaningless macro variable: 
		% adding up the micro one time using the raw values of X, one time using the exponential of X
		macro = zeros(1, npoints);
		for k = 1:(size(X,1));
			macro = macro + X(k,:);
		end 		
		
		% shuffled phases (shuffle each row)
		shuffled_micro = shuffle_rows(X);
		shuffled_macro = shuffle_rows(macro);
		
		% check distributions of a subset of parameters
		%{
		variables = {X, macro, shuffled_micro, shuffled_macro};
		titles = {['nodes'], ['sum of nodes'], ['shuffled nodes'], ['shuffled sum of nodes']};
		filenames = {['node'], ['sum'], ['shuffled_node'], ['shuffled_sum']};
		
		starting_point = npoints/2;
		end_point = (npoints/2)+100;
		
		A_str = param2str(coupling_vec(i));
		
		if (((i == 1) && (j == 1)) || ((i == 3) && (j == 3)) || ((i == 7) && (j == 7)) || ((i == 10) && (j == 10)));
			
			for h = 1:length(variables);
				variable = variables{h};
				
				num_var = size(variable, 1);
				r = randi([1 num_var],1,3);
				
				for k = 1:length(r);
					figure;
					histogram(variable(r(k),:)', nbins);
					title({['distribution of ', titles{h}]  ['var #', num2str(r(k)), ', A = ', num2str(coupling_vec(i)), ', beta = ', num2str(error_vec(j))]});
					ylabel('frequency');
					xlabel(['values of ', titles{h}, ', var #', num2str(r(k))]);
					exportgraphics(gcf, [pathout_plots network '_distr_' filenames{h} '_' num2str(r(k)) '_' A_str '_' num2str(npoints) '.png']);
				end
				
				for k = 1:length(r);
					figure;
					plot(variable(r(k),starting_point:end_point), 'LineWidth', 3);
					title({['time-series of ', titles{h}]  ['var #', num2str(r(k)), ', A = ', num2str(coupling_vec(i)), ', beta = ', num2str(error_vec(j))]});
					ylabel('value');
					xlabel('time');
					ylim([min(variable(r(k),starting_point:end_point))-0.3 max(variable(r(k),starting_point:end_point))+0.3]);
					xticklabels({num2str((npoints/2)) num2str(((npoints/2)+100)) num2str(((npoints/2)+200)) num2str(((npoints/2)+300)) num2str(((npoints/2)+400)) num2str(((npoints/2)+500))});
					exportgraphics(gcf, [pathout_plots network '_time_series_' filenames{h} '_' num2str(r(k)) '_' A_str '_' num2str(npoints) '.png']);
				end
				
				close all;
				
			end
		end 
		%}
		
		% { 
		synergy_capacity_practical_linear(i,j) = EmergencePsi(X', macro');
		downward_causation_practical_linear(i,j) = EmergenceDelta(X', macro');
		causal_decoupling_practical_linear(i,j) = synergy_capacity_practical_linear(i,j) - downward_causation_practical_linear(i,j);	
		
		synergy_capacity_practical_shuffled_micro(i,j) = EmergencePsi(shuffled_micro', macro');
		downward_causation_practical_shuffled_micro(i,j) = EmergenceDelta(shuffled_micro', macro');
		causal_decoupling_practical_shuffled_micro(i,j) = synergy_capacity_practical_shuffled_micro(i,j) - downward_causation_practical_shuffled_micro(i,j);

		synergy_capacity_practical_shuffled_macro(i,j) = EmergencePsi(X', shuffled_macro');
		downward_causation_practical_shuffled_macro(i,j) = EmergenceDelta(X', shuffled_macro');
		causal_decoupling_practical_shuffled_macro(i,j) = synergy_capacity_practical_shuffled_macro(i,j) - downward_causation_practical_shuffled_macro(i,j);

		synergy_capacity_practical_shuffled_micro_macro(i,j) = EmergencePsi(shuffled_micro', shuffled_macro');
		downward_causation_practical_shuffled_micro_macro(i,j) = EmergenceDelta(shuffled_micro', shuffled_macro');
		causal_decoupling_practical_shuffled_micro_macro(i,j) = synergy_capacity_practical_shuffled_micro_macro(i,j) - downward_causation_practical_shuffled_micro_macro(i,j);
		%}
		
		% dynamical dependence
		% store micro and macro variables in two structs (variables are transposed, so have time-points in rows)
		macro_variables.macro = macro';
		macro_variables.shuffled_macro = shuffled_macro';
		
		micro_variables.X = X';
		micro_variables.shuffled_micro = shuffled_micro';
		
		% get DD for all combinations of micro and top-level macro variables
		DD = get_DD(micro_variables, macro_variables, method_DD, tau, tau_steps);
		
		% extract DD for different combinations of micro and macro variables, and store it in arrays for
		% different parameter combinations of beta and A
		dd_X_macro(i,j) = DD.X_macro;
		dd_X_shuffled_macro(i,j) = DD.X_shuffled_macro;
		dd_shuffled_micro_macro(i,j) = DD.X_shuffled_macro;
		dd_shuffled_micro_shuffled_macro(i,j) = DD.shuffled_micro_shuffled_macro;
		
		clear DD;
		clear macro_variables;
		clear micro_variables;
		
		% average covariance/correlation matrix
		%{
		cov_X = cov(X');
		all_average_cov_X(i,j) = mean(nonzeros(tril(cov_X,-1)), 'all');

		corr_X = corrcov(cov_X);
		all_average_corr_X(i,j) = mean(nonzeros(tril(corr_X,-1)), 'all');
		%}
	end 
		
end 

%save([pathout_data network '_all_average_corr_X_' sim_index '.mat'], 'all_average_corr_X');
%save([pathout_data network '_all_average_cov_X_' sim_index '.mat'], 'all_average_cov_X');

%% calculating information atoms

%{
for i = 1:size(coupling_matrices, 3);
	
	coupling_matrix = coupling_matrices(:,:,i);
	disp(i)
	
	for j = 1:length(error_vec)

		err = error_vec(j);

		[X, Y] = sim_method(coupling_matrix, npoints, tau, err, corr);

		% PhiID 
		phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X(1,:), X(2,:), Y(1,:), Y(2,:), 'mmi'))';
		phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X(1,:), X(2,:), Y(1,:), Y(2,:), 'ccs'))';

	end 
		
end 

%}

%% storing information atoms & practical measures for different macro variables in struct files

% practical measures for different macro variables
% { 
emergence_practical = [];

emergence_practical.synergy_capacity_practical_linear = synergy_capacity_practical_linear;
emergence_practical.causal_decoupling_practical_linear = causal_decoupling_practical_linear;
emergence_practical.downward_causation_practical_linear = downward_causation_practical_linear;

emergence_practical.synergy_capacity_practical_shuffled_micro = synergy_capacity_practical_shuffled_micro;
emergence_practical.causal_decoupling_practical_shuffled_micro = causal_decoupling_practical_shuffled_micro;
emergence_practical.downward_causation_practical_shuffled_micro = downward_causation_practical_shuffled_micro;

emergence_practical.synergy_capacity_practical_shuffled_macro = synergy_capacity_practical_shuffled_macro;
emergence_practical.causal_decoupling_practical_shuffled_macro = causal_decoupling_practical_shuffled_macro;
emergence_practical.downward_causation_practical_shuffled_macro = downward_causation_practical_shuffled_macro;

emergence_practical.synergy_capacity_practical_shuffled_micro_macro = synergy_capacity_practical_shuffled_micro_macro;
emergence_practical.causal_decoupling_practical_shuffled_micro_macro = causal_decoupling_practical_shuffled_micro_macro;
emergence_practical.downward_causation_practical_shuffled_micro_macro = downward_causation_practical_shuffled_micro_macro;

save([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');

%}

% dynamical dependence
DD.dd_X_macro = dd_X_macro;
DD.dd_X_shuffled_macro = dd_X_shuffled_macro;
DD.dd_shuffled_micro_macro = dd_X_shuffled_macro;
DD.dd_shuffled_micro_shuffled_macro = dd_shuffled_micro_shuffled_macro;
				
save([pathout_data network '_DD' sim_index '.mat'], 'DD');

%% information atoms: extract single atoms such that rows are the couplings, and columns the errors

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

load([pathout_data network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
load([pathout_data network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');
load([pathout_data network '_emergence_practical' sim_index '.mat'], 'emergence_practical');
load([pathout_data network '_DD' sim_index '.mat'], 'DD');
load([pathout_data network '_all_atoms_err_coup_ccs' sim_index '.mat'], 'all_atoms_err_coup_ccs');
load([pathout_data network '_all_atoms_err_coup_mmi' sim_index '.mat'], 'all_atoms_err_coup_mmi');
load([pathout_data network '_all_average_corr_X' sim_index '.mat'], 'all_average_corr_X');
load([pathout_data network '_all_average_cov_X' sim_index '.mat'], 'all_average_cov_X');


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
	y_axis = {'0.0045', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '',  '', '0.09',  '', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', ... 
		'0.18', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '',  '', '0.27', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '', ... 
		'0.36', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '0.45', ''};
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

	xticks(1:100);
	yticks(1:100);
	set(gca,'TickLength',[0 0])
	yticklabels(y_axis);
	xticklabels(x_axis);
	
	if sim_index == '3';
		ylabel('zero coupling');
		xlabel('zero noise correlation');
	else 
		xlabel('noise correlation')%% calculating information atoms & practical measures
	end
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);
	ylabel('coupling strength');
end
%}

% synergistic capacity, downward causation, causal decoupling

% {
% heatmaps using imagesc()

% causal emergence
atoms = {emergence_ccs.synergy_capacity_ccs, emergence_mmi.synergy_capacity_mmi, ...
	emergence_practical.synergy_capacity_practical_linear, ...
	emergence_practical.synergy_capacity_practical_shuffled_micro, ...
	emergence_practical.synergy_capacity_practical_shuffled_macro, ...
	emergence_practical.synergy_capacity_practical_shuffled_micro_macro, ...
	emergence_ccs.downward_causation_ccs, emergence_mmi.downward_causation_mmi, ...
	emergence_practical.downward_causation_practical_linear, ...
	emergence_practical.downward_causation_practical_shuffled_micro, ...
	emergence_practical.downward_causation_practical_shuffled_macro, ...
	emergence_practical.downward_causation_practical_shuffled_micro_macro, ...
	emergence_ccs.causal_decoupling_ccs, emergence_mmi.causal_decoupling_mmi, ...
	emergence_practical.causal_decoupling_practical_linear, ...
	emergence_practical.causal_decoupling_practical_shuffled_micro, ...
	emergence_practical.causal_decoupling_practical_shuffled_macro, ...
	emergence_practical.causal_decoupling_practical_shuffled_micro_macro, ...
	all_average_cov_X, ...
	all_average_corr_X};

file_names = {'_all_err_coup_ccs_synergy_capacity', ...
	'_all_err_coup_mmi_synergy_capacity',  ...
	'_all_err_coup_synergy_capacity_practical_linear', ...
	'_all_err_coup_synergy_capacity_practical_shuffled_micro', ...
	'_all_err_coup_synergy_capacity_practical_shuffled_macro', ...
	'_all_err_coup_synergy_capacity_practical_shuffled_micro_macro', ...
	'_all_err_coup_ccs_downward_causation', ...
	'_all_err_coup_mmi_downward_causation', ...
	'_all_err_coup_downward_causation_practical_linear', ...
	'_all_err_coup_downward_causation_practical_shuffled_micro', ...
	'_all_err_coup_downward_causation_practical_shuffled_macro', ...
	'_all_err_coup_downward_causation_practical_shuffled_micro_macro', ...
	'_all_err_coup_ccs_causal_decoupling', ...
	'_all_err_coup_mmi_causal_decoupling', ...
	'_all_err_coup_causal_decoupling_practical_linear', ...
	'_all_err_coup_causal_decoupling_practical_shuffled_micro', ...
	'_all_err_coup_causal_decoupling_practical_shuffled_macro', ...
	'_all_err_coup_causal_decoupling_practical_shuffled_micro_macro', ...
	'_all_err_coup_average_cov_X', ...
	'_all_err_coup_average_corr_X'};

titles = {'synergy capacity ccs', 'synergy capacity mmi', ...
	'synergy capacity practical linear', ...
	'synergy capacity practical shuffled micro', ...
	'synergy capacity practical shuffled macro', ...
	'synergy capacity practical shuffled micro & macro', ...
	'downward causation ccs', 'downward causation mmi', ...
	'downward causation practical linear', ...
	'downward causation practical shuffled micro', ...
	'downward causation practical shuffled macro', ...
	'downward causation practical shuffled micro & macro', ...
	'causal decoupling ccs', ...
	'causal decoupling mmi', ...
	'causal decoupling practical linear', ...
	'causal decoupling practical shuffled micro', ...
	'causal decoupling practical shuffled macro', ...
	'causal decoupling practical shuffled micro & macro', ...
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
		ylabel('coupling strength');
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end
%}

close all;

% dynamical dependence

DD.dd_X_macro = dd_X_macro;
DD.dd_X_shuffled_macro = dd_X_shuffled_macro;
DD.dd_shuffled_micro_macro = dd_X_shuffled_macro;
DD.dd_shuffled_micro_shuffled_macro = dd_shuffled_micro_shuffled_macro;

atoms = {DD.dd_X_macro, ...
	DD.dd_X_shuffled_macro, ...
	DD.dd_shuffled_micro_macro, ...
	DD.dd_shuffled_micro_shuffled_macro};

file_names = {'_all_err_coup_dd_X_macro', ...
	'_all_err_coup_dd_X_shuffled_macro',  ...
	'_all_err_coup_dd_shuffled_micro_macro', ...
	'_all_err_coup_dd_shuffled_micro_shuffled_macro'};

titles = {'dynamical dependence micro & macro', ...
	'dynamical dependence micro & shuffled macro', ...
	'dynamical dependence shuffled micro & macro', ...
	'dynamical dependence shuffled micro & shuffled macro'};

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
		ylabel('coupling strength');
	end 
	
	title(titles{i});
	exportgraphics(gcf, [pathout_plots network file_names{i} sim_index '.png']);

end