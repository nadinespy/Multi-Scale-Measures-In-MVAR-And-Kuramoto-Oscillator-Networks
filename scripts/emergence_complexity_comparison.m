%% TO DO
% - where/how to store model information and parameter values for noise correlation & coupling strength in the saved mat-file?
% - fill struct file in a loop?

%%

clear all;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
% addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/PhiIDComparison/scripts/heatmaps'

PATHOUT1 = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/results/';
PATHOUT2 = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/results/plots/';
javaaddpath('infodynamics.jar');

%% choice of parameters

% time-lag and number of data points in time-series (same for all simulations)
npoints = 2000;
tau = 1;

% simulation method (options: statdata_corr_errors1(), statdata_corr_errors2())
sim_method = @statdata_corr_errors2;
% save plots & matrices according to simulation method (options: '1' (for statdata_corr_errors1()), '2' (for statdata_corr_errors2()))
sim_index = '2';

% network (options: '2node' for 2-node network with 100 different coupling strengths & noise correlations, '8node' for 8-node networks with different architectures)
network = '2node';
% error correlations (options: 100 different values for 2-node network, 6 different values for 8-node network)
error_vec = linspace(0.01, 0.9, 100); 

%% create coupling matrices depending on the chosen network size

if network == '2node'
	coupling_vec = linspace(0.01,0.45, 100);

	coupling_matrices = []; 
	for i = 1:length(coupling_vec)
		coupling_matrices(:,:,i) = coupling_vec(i)*ones(2);
	end 
	
else 
	load('nets.mat');
	net_names = fields(nets);
	nb_nets = length(net_names);

	coupling_matrices = [];
	for i = 1:size(net_names,1);
		coupling_matrices(:,:,i) = nets.(net_names{i});
	end 

end 

%% calculating information atoms

% instantiate variables to store atoms for different coupling matrices and noise correlations
phiid_all_err_coup_mmi = zeros(16, size(coupling_matrices,3), size(error_vec, 2));
phiid_all_err_coup_ccs = zeros(16, size(coupling_matrices,3), size(error_vec, 2));

for i = 1:size(coupling_matrices,3)
	
	coupling_matrix = coupling_matrices(:,:,i);
	spectral_radius = max(abs(eig(coupling_matrix)));
	
	for j = 1:length(error_vec)
		
		err = error_vec(j);
		X = sim_method(coupling_matrix, npoints, tau, err);
		phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
		phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
		phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		
	end 
	i
	%spectral_radius
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

save([PATHOUT1 network '_all_atoms_err_coup_ccs' sim_index '.mat'],'all_atoms_err_coup_ccs');
save([PATHOUT1 network '_all_atoms_err_coup_mmi' sim_index '.mat'],'all_atoms_err_coup_mmi');

%% double-redundancy & double-synergy

% heatmaps

if network == '2node'
	x_axis = {'0.09', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''};
	y_axis = {'0.0045', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '',  '', '0.09',  '', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', ... 
		'0.18', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '',  '', '0.27', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '', ... 
		'0.36', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '0.45', ''};
else 
	x_axis = {'0.15', '0.3', '0.45', '0.6', '0.75', '0.9'};
	y_axis = {'fully connected', 'optimal A', 'optimal B', 'ring', 'small world', 'uni ring'};
end 

% using matlab built-in function for heatmaps

atoms = {all_atoms_err_coup_ccs.rtr, all_atoms_err_coup_mmi.rtr, all_atoms_err_coup_ccs.sts, all_atoms_err_coup_mmi.sts};
file_names = {'_all_err_coup_ccs_rtr', '_all_err_coup_mmi_rtr', '_all_err_coup_ccs_sts', '_all_err_coup_mmi_sts'};
titles = {'double-redundancy ccs', 'double-redundancy mmi', 'double-synergy ccs', 'double-synergy mmi'};

for i = 1:size(atoms,2)
	
	figure;
	clf
	h = heatmap(atoms{i}, 'Colormap', parula, 'ColorbarVisible', 'on', 'CellLabelColor', 'none') ;
	h.YDisplayLabels = y_axis;
	h.XDisplayLabels = x_axis;
	h.XLabel = 'noise correlation';
	
	if network == '2node'
		h.YLabel = 'coupling strength';
	else 
		h.YLabel = 'network architecture';
	end
	
	title(titles{i});
	exportgraphics(gcf, [PATHOUT2 network file_names{i} sim_index '.png']);

end


%% synergistic/emergent capacity, downward causation, causal decoupling

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
emergence_css.causal_decoupling_ccs = causal_decoupling_ccs;
emergence_css.downward_causation_ccs = downward_causation_ccs;

emergence_mmi.synergy_capacity_mmi = synergy_capacity_mmi;
emergence_mmi.causal_decoupling_mmi = causal_decoupling_mmi;
emergence_mmi.downward_causation_mmi = downward_causation_mmi;

save([PATHOUT1 network '_emergence_ccs' sim_index '.mat'], 'emergence_ccs');
save([PATHOUT1 network '_emergence_mmi' sim_index '.mat'], 'emergence_mmi');

% heatmaps using matlab built-in function for heatmaps:

atoms = {synergy_capacity_ccs, synergy_capacity_mmi, downward_causation_ccs, downward_causation_mmi, causal_decoupling_ccs, causal_decoupling_mmi};
file_names = {'_all_err_coup_ccs_synergy_capacity', '_all_err_coup_mmi_synergy_capacity', '_all_err_coup_ccs_downward_causation', '_all_err_coup_mmi_downward_causation', '_all_err_coup_ccs_causal_decoupling', '_all_err_coup_mmi_causal_decoupling'};
titles = {'synergy capacity ccs', 'synergy capacity mmi', 'downward causation ccs', 'downward causation mmi', 'causal decoupling ccs', 'causal decoupling mmi'};

for i = 1:size(atoms,2)
	
	figure;
	clf
	h = heatmap(atoms{i}, 'Colormap', parula, 'ColorbarVisible', 'on', 'CellLabelColor', 'none') ;
	h.YDisplayLabels = y_axis;
	h.XDisplayLabels = x_axis;
	h.XLabel = 'noise correlation';
	
	if network == '2node'
		h.YLabel = 'coupling strength';
	else 
		h.YLabel = 'network architecture';
	end
	
	title(titles{i});
	exportgraphics(gcf, [PATHOUT2 network file_names{i} sim_index '.png']);

end

close all;
