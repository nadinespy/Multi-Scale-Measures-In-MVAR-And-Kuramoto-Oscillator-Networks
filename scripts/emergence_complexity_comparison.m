clear all;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison'
addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/scripts'
addpath '/media/nadinespy/NewVolume/my_stuff/work/toolboxes_matlab'
% addpath '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/PhiIDComparison/scripts/heatmaps'

PATHOUT = '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/results/';
javaaddpath('infodynamics.jar');

%% calculating information atoms

% System size, coupling matrix, and vector of noise values
nvar = 2;
npoints = 2000;
tau = 1;
error_vec   = linspace(0.01, 0.99, 100);
coupling_vec = linspace(0.01,0.45, 100);

phiid_all_err_coup_mmi = zeros(16, size(coupling_vec,2), size(error_vec, 2));
phiid_all_err_coup_ccs = zeros(16, size(coupling_vec,2), size(error_vec, 2));


for i = 1:length(coupling_vec)
	A = coupling_vec(i)*ones(nvar);
	
	for j = 1:length(error_vec)
		err = error_vec(j);
		spectral_radius = max(abs(eig(A)));
		%X = statdata_corr_errors(A, npoints, err);
		X = statdata_corr_errors2(A, npoints, err);
		%X = statdata_corr_errors3(A, npoints, err);
		phiid_all_err_coup_mmi(:,i,j) = struct2array(PhiIDFull(X, tau, 'MMI'))';
		phiid_all_err_coup_ccs(:,i,j) = struct2array(PhiIDFull(X, tau, 'ccs'))';
		
	end 
	i
	%spectral_radius
end 

save([PATHOUT 'phiid_all_err_coup_ccs2.mat'],'phiid_all_err_coup_ccs');
save([PATHOUT 'phiid_all_err_coup_mmi2.mat'],'phiid_all_err_coup_mmi');

% save([PATHOUT 'phiid_all_err_coup_ccs_rtr.mat'],'phiid_all_err_coup_ccs_rtr');
% save([PATHOUT 'phiid_all_err_coup_ccs_sts.mat'],'phiid_all_err_coup_ccs_sts');
% save([PATHOUT 'phiid_all_err_coup_mmi_rtr.mat'],'phiid_all_err_coup_ccs_mmi');
% save([PATHOUT 'phiid_all_err_coup_mmi_sts.mat'],'phiid_all_err_coup_ccs_mmi');

%% double-redundancy & double-synergy

% extract double-redundancy & synergy term
% rows: couplings; errors: columns
phiid_all_err_coup_ccs_rtr = squeeze(phiid_all_err_coup_ccs(1,:,:));			% coupling vector will now be in the rows, and error vector in the columns
phiid_all_err_coup_ccs_sts = squeeze(phiid_all_err_coup_ccs(16,:,:));

phiid_all_err_coup_mmi_rtr = squeeze(phiid_all_err_coup_mmi(1,:,:));
phiid_all_err_coup_mmi_sts = squeeze(phiid_all_err_coup_mmi(16,:,:));

% heatmaps
x_axis = {'0.01', '', '', '', '', '', '', '', '', '0.1', '', '', '', '', '', '', '', '', '', '0.2', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', '', '', '', ... 
	'0.4', '', '', '', '', '', '', '', '', '', '0.5', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '0.7', '', '', '', '', '', '', '', '', '', ... 
	'0.8', '', '', '', '', '', '', '', '', '', '0.9', '', '', '', '', '', '', '', '', '0.99', ''}; 
y_axis = {'0.01', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '',  '', '0.1',  '', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', ... 
	'0.2', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '',  '', '0.3', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '', ... 
	'0.4', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '',  '0.5'};

% using matlab built-in function for heatmaps:
figure;
clf
h = heatmap(phiid_all_err_coup_mmi_rtr, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('double-redundancy mmi');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_mmi_rtr2.png']);

figure;
clf
h = heatmap(phiid_all_err_coup_mmi_sts, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('double-synergy mmi');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_mmi_sts2.png']);

figure;
clf
h = heatmap(phiid_all_err_coup_ccs_rtr, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('double-redundancy ccs');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_ccs_rtr2.png']);

figure;
clf
h = heatmap(phiid_all_err_coup_ccs_sts, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('double-synergy ccs');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_ccs_sts2.png']);

close all;

%% synergistic capacity, downward causation, causal decoupling

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
phiid_all_err_coup_ccs_rtx = squeeze(phiid_all_err_coup_ccs(2,:,:));			% coupling vector will now be in the rows, and error vector in the columns
phiid_all_err_coup_ccs_rty = squeeze(phiid_all_err_coup_ccs(3,:,:));
phiid_all_err_coup_ccs_rts = squeeze(phiid_all_err_coup_ccs(4,:,:));
phiid_all_err_coup_ccs_xtr = squeeze(phiid_all_err_coup_ccs(5,:,:));
phiid_all_err_coup_ccs_xtx = squeeze(phiid_all_err_coup_ccs(6,:,:));
phiid_all_err_coup_ccs_xty = squeeze(phiid_all_err_coup_ccs(7,:,:));
phiid_all_err_coup_ccs_xts = squeeze(phiid_all_err_coup_ccs(8,:,:));
phiid_all_err_coup_ccs_ytr = squeeze(phiid_all_err_coup_ccs(9,:,:));
phiid_all_err_coup_ccs_ytx = squeeze(phiid_all_err_coup_ccs(10,:,:));
phiid_all_err_coup_ccs_yty = squeeze(phiid_all_err_coup_ccs(11,:,:));
phiid_all_err_coup_ccs_yts = squeeze(phiid_all_err_coup_ccs(12,:,:));
phiid_all_err_coup_ccs_str = squeeze(phiid_all_err_coup_ccs(13,:,:));
phiid_all_err_coup_ccs_stx = squeeze(phiid_all_err_coup_ccs(14,:,:));
phiid_all_err_coup_ccs_sty = squeeze(phiid_all_err_coup_ccs(15,:,:));

phiid_all_err_coup_mmi_rtx = squeeze(phiid_all_err_coup_mmi(2,:,:));			% coupling vector will now be in the rows, and error vector in the columns
phiid_all_err_coup_mmi_rty = squeeze(phiid_all_err_coup_mmi(3,:,:));
phiid_all_err_coup_mmi_rts = squeeze(phiid_all_err_coup_mmi(4,:,:));
phiid_all_err_coup_mmi_xtr = squeeze(phiid_all_err_coup_mmi(5,:,:));
phiid_all_err_coup_mmi_xtx = squeeze(phiid_all_err_coup_mmi(6,:,:));
phiid_all_err_coup_mmi_xty = squeeze(phiid_all_err_coup_mmi(7,:,:));
phiid_all_err_coup_mmi_xts = squeeze(phiid_all_err_coup_mmi(8,:,:));
phiid_all_err_coup_mmi_ytr = squeeze(phiid_all_err_coup_mmi(9,:,:));
phiid_all_err_coup_mmi_ytx = squeeze(phiid_all_err_coup_mmi(10,:,:));
phiid_all_err_coup_mmi_yty = squeeze(phiid_all_err_coup_mmi(11,:,:));
phiid_all_err_coup_mmi_yts = squeeze(phiid_all_err_coup_mmi(12,:,:));
phiid_all_err_coup_mmi_str = squeeze(phiid_all_err_coup_mmi(13,:,:));
phiid_all_err_coup_mmi_stx = squeeze(phiid_all_err_coup_mmi(14,:,:));
phiid_all_err_coup_mmi_sty = squeeze(phiid_all_err_coup_mmi(15,:,:));

% synergy
synergy_capacity_ccs = phiid_all_err_coup_ccs_str + ...
	phiid_all_err_coup_ccs_stx + phiid_all_err_coup_ccs_sty + phiid_all_err_coup_ccs_sts;

downward_causation_ccs = phiid_all_err_coup_ccs_str + phiid_all_err_coup_ccs_stx + phiid_all_err_coup_ccs_sty;

synergy_capacity_mmi = phiid_all_err_coup_mmi_str + ...
	phiid_all_err_coup_mmi_stx + phiid_all_err_coup_mmi_sty + phiid_all_err_coup_mmi_sts;

downward_causation_mmi = phiid_all_err_coup_mmi_str + phiid_all_err_coup_mmi_stx + phiid_all_err_coup_mmi_sty;

causal_decoupling_ccs = synergy_capacity_ccs - downward_causation_ccs;
causal_decoupling_mmi = synergy_capacity_mmi - downward_causation_mmi;

% heatmaps
x_axis = {'0.01', '', '', '', '', '', '', '', '', '0.1', '', '', '', '', '', '', '', '', '', '0.2', '', '', '', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', '', '', '', ... 
	'0.4', '', '', '', '', '', '', '', '', '', '0.5', '', '', '', '', '', '', '', '', '', '0.6', '', '', '', '', '', '', '', '', '', '0.7', '', '', '', '', '', '', '', '', '', ... 
	'0.8', '', '', '', '', '', '', '', '', '', '0.9', '', '', '', '', '', '', '', '', '0.99', ''}; 
y_axis = {'0.01', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '',  '', '0.1',  '', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', ... 
	'0.2', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '',  '', '0.3', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '', ... 
	'0.4', '', '', '', '', '', '', '', '', '',  '', '', '', '', '', '', '', '', '', '',  '0.5'};

% using matlab built-in function for heatmaps:
figure;
clf
h = heatmap(synergy_capacity_ccs, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('synergy capacity ccs');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_ccs_synergy_capacity2.png']);

figure;
clf
h = heatmap(synergy_capacity_mmi, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('synergy capacity mmi');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_mmi_synergy_capacity2.png']);

figure;
clf
h = heatmap(downward_causation_ccs, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('downward causation ccs');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_ccs_downward_causation2.png']);

figure;
clf
h = heatmap(downward_causation_mmi, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('downward causation mmi');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_mmi_downward_causation2.png']);

figure;
clf
h = heatmap(causal_decoupling_mmi, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('causal decoupling mmi');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_mmi_causal_decoupling2.png']);

figure;
clf
h = heatmap(causal_decoupling_ccs, 'Colormap', parula, 'ColorbarVisible', 'on') ;
h.YDisplayLabels = y_axis;
h.XDisplayLabels = x_axis;
h.XLabel = 'noise correlation';
h.YLabel = 'coupling strength';
title('causal decoupling ccs');
exportgraphics(gcf, [PATHOUT '2node_phiid_all_err_coup_ccs_causal_decoupling2.png']);

close all;
