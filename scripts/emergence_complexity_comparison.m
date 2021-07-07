clear all;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
% addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/PhiIDComparison/scripts/heatmaps'

PATHOUT = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/results/';
javaaddpath('infodynamics.jar');

% choose simulation method (options: statdata_corr_errors1(), statdata_corr_errors2())
sim_method = @statdata_corr_errors1;

% save plots & matrices according to simulation method (options: '1' (for statdata_corr_errors1()), '2' (for statdata_corr_errors2()))
sim_index = '1';


%% calculating information atoms

% System size, coupling matrix, time-lag, and vector of noise values
nvar = 2;
npoints = 2000;
tau = 1;
error_vec   = linspace(0.01, 0.9, 100);
coupling_vec = linspace(0.01,0.45, 100);

phiid_all_err_coup_mmi = zeros(16, size(coupling_vec,2), size(error_vec, 2));
phiid_all_err_coup_ccs = zeros(16, size(coupling_vec,2), size(error_vec, 2));


for i = 1:length(coupling_vec)
	A = coupling_vec(i)*ones(nvar);
	
	for j = 1:length(error_vec)
		err = error_vec(j);
		spectral_radius = max(abs(eig(A)));
		%X = statdata_corr_errors(A, npoints, tau, err);
		X = sim_method(A, npoints, tau, err);
		phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
		phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		
	end 
	i
	%spectral_radius
end 

% allocating variable names for the atoms in a struct
all_atoms_err_coup_mmi = [];
all_atoms_err_coup_mmi.rtr = phiid_all_err_coup_mmi(1,:,:);
all_atoms_err_coup_mmi.rtx = phiid_all_err_coup_mmi(2,:,:);
all_atoms_err_coup_mmi.rty = phiid_all_err_coup_mmi(3,:,:);
all_atoms_err_coup_mmi.rts = phiid_all_err_coup_mmi(4,:,:);
all_atoms_err_coup_mmi.xtr = phiid_all_err_coup_mmi(5,:,:);
all_atoms_err_coup_mmi.xtx = phiid_all_err_coup_mmi(6,:,:);
all_atoms_err_coup_mmi.xty = phiid_all_err_coup_mmi(7,:,:);
all_atoms_err_coup_mmi.xts = phiid_all_err_coup_mmi(8,:,:);
all_atoms_err_coup_mmi.ytr = phiid_all_err_coup_mmi(9,:,:);
all_atoms_err_coup_mmi.ytx = phiid_all_err_coup_mmi(10,:,:);
all_atoms_err_coup_mmi.yty = phiid_all_err_coup_mmi(11,:,:);
all_atoms_err_coup_mmi.yts = phiid_all_err_coup_mmi(12,:,:);
all_atoms_err_coup_mmi.str = phiid_all_err_coup_mmi(13,:,:);
all_atoms_err_coup_mmi.stx = phiid_all_err_coup_mmi(14,:,:);
all_atoms_err_coup_mmi.sty = phiid_all_err_coup_mmi(15,:,:);
all_atoms_err_coup_mmi.sts = phiid_all_err_coup_mmi(16,:,:);

all_atoms_err_coup_ccs = [];
all_atoms_err_coup_ccs.rtr = phiid_all_err_coup_ccs(1,:,:);
all_atoms_err_coup_ccs.rtx = phiid_all_err_coup_ccs(2,:,:);
all_atoms_err_coup_ccs.rty = phiid_all_err_coup_ccs(3,:,:);
all_atoms_err_coup_ccs.rts = phiid_all_err_coup_ccs(4,:,:);
all_atoms_err_coup_ccs.xtr = phiid_all_err_coup_ccs(5,:,:);
all_atoms_err_coup_ccs.xtx = phiid_all_err_coup_ccs(6,:,:);
all_atoms_err_coup_ccs.xty = phiid_all_err_coup_ccs(7,:,:);
all_atoms_err_coup_ccs.xts = phiid_all_err_coup_ccs(8,:,:);
all_atoms_err_coup_ccs.ytr = phiid_all_err_coup_ccs(9,:,:);
all_atoms_err_coup_ccs.ytx = phiid_all_err_coup_ccs(10,:,:);
all_atoms_err_coup_ccs.yty = phiid_all_err_coup_ccs(11,:,:);
all_atoms_err_coup_ccs.yts = phiid_all_err_coup_ccs(12,:,:);
all_atoms_err_coup_ccs.str = phiid_all_err_coup_ccs(13,:,:);
all_atoms_err_coup_ccs.stx = phiid_all_err_coup_ccs(14,:,:);
all_atoms_err_coup_ccs.sty = phiid_all_err_coup_ccs(15,:,:);
all_atoms_err_coup_ccs.sts = phiid_all_err_coup_ccs(16,:,:);

save([PATHOUT 'all_atoms_err_coup_ccs' sim_index '.mat'],'all_atoms_err_coup_ccs');
save([PATHOUT 'all_atoms_err_coup_mmi' sim_index '.mat'],'all_atoms_err_coup_mmi');

%% double-redundancy & double-synergy

% extract single atoms (double-redundancy & synergy term)
% such that rows are the couplings, and columns the errors
all_err_coup_ccs_rtr = squeeze(all_atoms_err_coup_ccs.rtr(1,:,:));			% coupling vector will now be in the rows, and error vector in the columns
all_err_coup_ccs_sts = squeeze(all_atoms_err_coup_ccs.sts(1,:,:));

all_err_coup_mmi_rtr = squeeze(all_atoms_err_coup_mmi.rtr(1,:,:));			% coupling vector will now be in the rows, and error vector in the columns
all_err_coup_mmi_sts = squeeze(all_atoms_err_coup_mmi.sts(1,:,:));

% heatmaps
x_axis = {'0.09', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', ... 
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', '', ... 
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '0.9', ''}; 
y_axis = {'0.0045', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '',  '', '0.09',  '', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', ... 
	'0.18', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '',  '', '0.27', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '', ... 
	'0.36', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '0.45', ''};

% using matlab built-in function for heatmaps:

atoms = {all_err_coup_ccs_rtr, all_err_coup_mmi_rtr, all_err_coup_ccs_sts, all_err_coup_mmi_sts};
file_names = {'all_err_coup_ccs_rtr', 'all_err_coup_mmi_rtr', 'all_err_coup_ccs_sts', 'all_err_coup_mmi_sts'};
titles = {'double-redundancy ccs', 'double-redundancy mmi', 'double-synergy ccs', 'double-synergy mmi'};

for i = 1:size(atoms,2)
	
	figure;
	clf
	h = heatmap(atoms{i}, 'Colormap', parula, 'ColorbarVisible', 'on') ;
	h.YDisplayLabels = y_axis;
	h.XDisplayLabels = x_axis;
	h.XLabel = 'noise correlation';
	h.YLabel = 'coupling strength';
	title(titles{i});
	exportgraphics(gcf, [PATHOUT '2node_' file_names{i} sim_index '.png']);

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

% rows: couplings; errors: columns
all_err_coup_ccs_rtx = squeeze(all_atoms_err_coup_ccs.rtx(1,:,:));			% coupling vector will now be in the rows, and error vector in the columns
all_err_coup_ccs_rty = squeeze(all_atoms_err_coup_ccs.rty(1,:,:));
all_err_coup_ccs_rts = squeeze(all_atoms_err_coup_ccs.rts(1,:,:));
all_err_coup_ccs_xtr = squeeze(all_atoms_err_coup_ccs.xtr(1,:,:));
all_err_coup_ccs_xtx = squeeze(all_atoms_err_coup_ccs.xtx(1,:,:));
all_err_coup_ccs_xty = squeeze(all_atoms_err_coup_ccs.xty(1,:,:));
all_err_coup_ccs_xts = squeeze(all_atoms_err_coup_ccs.xts(1,:,:));
all_err_coup_ccs_ytr = squeeze(all_atoms_err_coup_ccs.ytr(1,:,:));
all_err_coup_ccs_ytx = squeeze(all_atoms_err_coup_ccs.ytx(1,:,:));
all_err_coup_ccs_yty = squeeze(all_atoms_err_coup_ccs.yty(1,:,:));
all_err_coup_ccs_yts = squeeze(all_atoms_err_coup_ccs.yts(1,:,:));
all_err_coup_ccs_str = squeeze(all_atoms_err_coup_ccs.str(1,:,:));
all_err_coup_ccs_stx = squeeze(all_atoms_err_coup_ccs.stx(1,:,:));
all_err_coup_ccs_sty = squeeze(all_atoms_err_coup_ccs.sty(1,:,:));

all_err_coup_mmi_rtx = squeeze(all_atoms_err_coup_mmi.rtx(1,:,:));			% coupling vector will now be in the rows, and error vector in the columns
all_err_coup_mmi_rty = squeeze(all_atoms_err_coup_mmi.rty(1,:,:));
all_err_coup_mmi_rts = squeeze(all_atoms_err_coup_mmi.rts(1,:,:));
all_err_coup_mmi_xtr = squeeze(all_atoms_err_coup_mmi.xtr(1,:,:));
all_err_coup_mmi_xtx = squeeze(all_atoms_err_coup_mmi.xtx(1,:,:));
all_err_coup_mmi_xty = squeeze(all_atoms_err_coup_mmi.xty(1,:,:));
all_err_coup_mmi_xts = squeeze(all_atoms_err_coup_mmi.xts(1,:,:));
all_err_coup_mmi_ytr = squeeze(all_atoms_err_coup_mmi.ytr(1,:,:));
all_err_coup_mmi_ytx = squeeze(all_atoms_err_coup_mmi.ytx(1,:,:));
all_err_coup_mmi_yty = squeeze(all_atoms_err_coup_mmi.yty(1,:,:));
all_err_coup_mmi_yts = squeeze(all_atoms_err_coup_mmi.yts(1,:,:));
all_err_coup_mmi_str = squeeze(all_atoms_err_coup_mmi.str(1,:,:));
all_err_coup_mmi_stx = squeeze(all_atoms_err_coup_mmi.stx(1,:,:));
all_err_coup_mmi_sty = squeeze(all_atoms_err_coup_mmi.sty(1,:,:));

% synergy
synergy_capacity_ccs = all_err_coup_ccs_str + ...
	all_err_coup_ccs_stx + all_err_coup_ccs_sty + all_err_coup_ccs_sts;

downward_causation_ccs = all_err_coup_ccs_str + all_err_coup_ccs_stx + all_err_coup_ccs_sty;

synergy_capacity_mmi = all_err_coup_mmi_str + ...
	all_err_coup_mmi_stx + all_err_coup_mmi_sty + all_err_coup_mmi_sts;

downward_causation_mmi = all_err_coup_mmi_str + all_err_coup_mmi_stx + all_err_coup_mmi_sty;

causal_decoupling_ccs = synergy_capacity_ccs - downward_causation_ccs;
causal_decoupling_mmi = synergy_capacity_mmi - downward_causation_mmi;

% heatmaps using matlab built-in function for heatmaps:

atoms = {synergy_capacity_ccs, synergy_capacity_mmi, downward_causation_ccs, downward_causation_mmi, causal_decoupling_ccs, causal_decoupling_mmi};
file_names = {'all_err_coup_ccs_synergy_capacity', 'all_err_coup_mmi_synergy_capacity', 'all_err_coup_ccs_downward_causation', 'all_err_coup_mmi_downward_causation', 'all_err_coup_ccs_causal_decoupling', 'all_err_coup_mmi_causal_decoupling'};
titles = {'synergy capacity ccs', 'synergy capacity mmi', 'downward causation ccs', 'downward causation mmi', 'causal decoupling ccs', 'causal decoupling mmi'};

for i = 1:size(atoms,2)
	
	figure;
	clf
	h = heatmap(atoms{i}, 'Colormap', parula, 'ColorbarVisible', 'on') ;
	h.YDisplayLabels = y_axis;
	h.XDisplayLabels = x_axis;
	h.XLabel = 'noise correlation';
	h.YLabel = 'coupling strength';
	title(titles{i});
	exportgraphics(gcf, [PATHOUT '2node_' file_names{i} sim_index '.png']);

end

close all;
