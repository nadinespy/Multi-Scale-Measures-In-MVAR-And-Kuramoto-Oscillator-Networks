clear all;
close all;
clc;

mvgc_path  = getenv('MVGC2_PATH');
run(fullfile(mvgc_path,'startup'));

cd '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/scripts'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/scripts/analytical/'
pathout_plots = {'/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/results/plots/2node_mvar/integration_analytical/'};

n = 2;		% number of variables in network
time_lag = 1;	% number of time-lags

%% CREATE CONNECTIVITY MATRIX

%{
load('all_nets.mat');
eight_node_net_names = fields(all_nets);
for i = 1:size(eight_node_net_names,1);
	net = all_nets.(eight_node_net_names{i});
	coupling_matrices(:,:,i) = net;
end
%}

%{										
% two-node matrices of the form [a c; 0 b] for a = 0.8; b = 0.7; 
% and c = y where y is one value from [-2:1000:2]
a = 0.8;% 	% positive noise correlation
% 	rho = sqrt(1-exp(-all_rmi(j)));
% 	R = [1 rho; rho 1];

b = 0.7;
cr = 2;
c = linspace(-cr,cr,1001)';
for i = 1:length(c);
	coupling_matrix = [a c(i); 0, b];
	coupling_matrices(:,:,i) = coupling_matrix;
end
%}

% {
% connectivity matrices with same couplings (including self-connection)
two_node_couplings	= linspace(0.0045,0.45,100);	% array of coupling values 
csfac				= linspace(0,0.99,100)';	% array of connectivity scale factors

% create two-node matrices of the form [a a; a a] for coupling value a
for i = 1:length(two_node_couplings);
	coupling_matrices(:,:,i) = two_node_couplings(i)*ones(2);
end  
%}

%%

% variation of residuals' mutual information and noise correlation
all_rmi			= linspace(0.0,1,100);		% array of rmi values 
noise_corrs			= linspace(0.0,0.9,100);	% array of noise correlations 
										
% selection of c-factors and couplings that are kept fixed for 
% different calculations/plots 
csfac_indices = [1 50 100];
coup_indices  = [1 50 100];

% choose whether correlations between noise sources should all 
% be positive ('positive'), negative ('negative'), or alternating 
% between positive and negative ('mixed')
signs_noise_corrs = 'positive'; 

% select specific rmi value for line plot
rmi_index = 10;

% calculate integration across values for rmi, noise correlations, 
% and couplings, fix c-factor
for k = 1:length(csfac_indices)
	csfac_index = csfac_indices(k);	% fix c-factor
	rmi_coup_fixed_csfac;			% calculate integration across rmi & couplings
	plot_rmi_coup_fixed_csfac;		% make plots
	clear I;
	noise_corr_coup_fixed_csfac;		% calculate integration across noise corrs & couplings
	plot_noise_corr_coup_fixed_csfac;	% make plots
	clear I;
end 

% calculate integration across values for rmi, noise correlations, 
% and c-factors, fix coupling
for k = 1:length(coup_indices)
	coup_index = coup_indices(k);		% fix coupling
	rmi_csfac_fixed_coup;			% calculate integration across rmi & c-factors 
	plot_rmi_csfac_fixed_coup;		% make plots
	clear I;
	clear phiidCE_MMI;
	clear phiidCE_CCS;
	clear phiidRed_MMI;
	clear phiidRed_CCS;
	clear phiidSyn_MMI;
	clear phiidSyn_CCS;
	noise_corr_csfac_fixed_coup;		% calculate integration across noise corrs & c-factors
	plot_noise_corr_csfac_fixed_coup;	% make plots
	clear I;
	clear phiidCE_MMI;
	clear phiidCE_CCS;
	clear phiidRed_MMI;
	clear phiidRed_CCS;
	clear phiidSyn_MMI;
	clear phiidSyn_CCS;
end 


