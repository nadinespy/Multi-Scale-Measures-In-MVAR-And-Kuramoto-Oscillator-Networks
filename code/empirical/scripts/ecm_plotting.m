%% PLOTTING FOR KURAMOTO OSCILLATORS WITH DIFFERENT COUPLINGS, BETAS, MACRO & MICRO VARIABLES

% TODO
% - make axes go from 0 at lower left corner (both x- and y-axis)

% SCRIPT DESCRIPTION

% major sections in this script:
%	- distributions of micro and macro variables						--> get_all_kuramoto_distr_plots()	
%
%	- heatmaps for correlations between micro & macro variables				--> get_all_kuramoto_corr_heatmaps();
%														    uses plot_heatmaps()
%
%	- heatmaps for emergence measures								--> get_all_emergence_heatmaps();
%
%	- scatter plots for practical CE, with fixed A, and varying beta			--> get_all_kuramoto_practCE_scatterplots();
%														    uses plot_scatterplos_measure()
%
%	- scatter plots for dynamical dependence, with fixed A, and varying beta	--> get_all_kuramoto_DD_scatterplots();
%														    uses plot_scatterplos_measure()
%		
%	- scatter plots for sigma met mean & sigma chi mean,					--> get_all_kuramoto_met_chi_scatterplots();
%       with fixed A, and varying beta								    uses plot_scatterplot_met_chi()

%% plotting parameter specification

% number of bins for plotting histograms of micro and macro variables
nbins = 100;								

%**************************************************************************
% KURAMOTO OSCILLATORS
%**************************************************************************

%{
% axes ticks for heatmaps with beta on x-axis, and A on y-axis
%x_axis_heatmaps = {'0.04', '', '', '0.16', '', '', '0.28', '', '', '0.4'};
x_axis_heatmaps = {'0.08', '', '', '0.32', '', '', '0.56', '', '', '0.8'};
y_axis_heatmaps = {'0.08', '', '', '0.32', '', '', '0.56', '', '', '0.8'};

y_label_heatmaps = 'A';
x_label_heatmaps = 'beta';

% scatterplots

% axes ticks for scatterplots with beta on x-axis, and Shannon-CE/DC/CD 
% on y-axis
y_label_practCE_scatterplots = 'Shannon-CE';
x_label_practCE_scatterplots = 'beta';

% axes ticks for scatterplots with beta on x-axis, and dynamical dependence 
% on y-axis
y_label_DD_scatterplots = 'dynamical dependence';
x_label_DD_scatterplots = 'beta';

% axes ticks for scatterplots with beta on x-axis, and PhiID-CE/DC/CD 
% on y-axis
y_label_DD_scatterplots = 'PhiID-CE';
x_label_DD_scatterplots = 'beta';

% axes ticks for scatterplots with beta on x-axis, and sigma chi/sigma met 
% on y-axis
y_label_chi_met = ' ';
x_label_chi_met = 'beta';
%}

%**************************************************************************
% MULTIVARIATE AUTOREGRESSIVE NETWORK 
%**************************************************************************

% axes ticks for heatmaps with beta on x-axis, and A on y-axis
%x_axis_heatmaps = {'0.04', '', '', '0.16', '', '', '0.28', '', '', '0.4'};
x_axis_heatmaps = {'0.09', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.32', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.56', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.9'};

y_axis_heatmaps = {'0.0045', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '0.15', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '0.3', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
	'', '', '', '', '', '', '', '', '', '0.45'};

y_label_heatmaps = 'coupling';
x_label_heatmaps = 'noise correlation';

%% load struct
struct_prefix = '';

% load struct...


%% plotting
%% distribution plots of micro and macro variables (selected nodes)

% to-do: make random selection of which nodes to plot

% loop over all values of npoints

%{
  
get_all_kuramoto_distr_plots(data, nbins, network, A_vec, beta_vec, all_npoints, ...
		pathout_data_sim_time_series, pathout_plots_distributions);
		
%}

%% heatmaps for metastability & chimera index

%{

get_km_met_chi_heatmaps(network, npoints, x_axis_heatmaps, y_axis_heatmaps, ...
	x_label_heatmaps, y_label_heatmaps, pathout_data_sync, pathout_plots_sigma_chi, ...
	pathout_plots_sigma_met)

%}
	
%% correlations between micro & macro variables

% loop over all values of npoints

%{

get_all_kuramoto_corr_heatmaps(network, all_npoints, x_axis_heatmaps, y_axis_heatmaps, ...
	x_label_heatmaps, y_label_heatmaps, pathout_data_mean_corr, pathout_plots_mean_corr)

%}
	
%% heatmaps for emergence measures

% make selection: choose all instances of, e. g., Shannon-based CE/DC/CD

%{
emergence_results_shannonCE_DC_CD = [];
for g = 1:length(emergence_results)
	if strfind(emergence_results(1,g).measure, 'shannon') == 1
		emergence_results_shannonCE_DC_CD = [emergence_results_shannonCE_DC_CD, emergence_results(1,g)];
	end 
end
%}

% make plots that have the same axes ticks and labels
get_all_emergence_heatmaps(network, emergence_results, x_axis_heatmaps, y_axis_heatmaps, ...
	x_label_heatmaps, y_label_heatmaps, pathout_plots_measures);

%}

	
%% scatter plots for practical CE for discretized variables, with fixed A, and varying beta

% loop over all values of npoints, tau

%{

get_all_kuramoto_practCE_scatterplots(network, A_vec, beta_vec, taus, all_npoints, x_label_practCE_scatterplots, ...
		y_label_practCE_scatterplots, pathout_data_pract_ce, pathout_plots_pract_ce_sigma_chi, ...
		pathout_plots_pract_ce_mean_pair_sync, pathout_plots_pract_ce_control_vars)
	
%} 

%% scatter plots for dynamical dependence with non-Gaussian, continuous variables, with fixed A, and varying beta

% loop over all values of npoints, tau

%{

get_all_kuramoto_DD_scatterplots(network, A_vec, beta_vec, taus, all_npoints, x_label_DD_scatterplots, ...
		y_label_DD_scatterplots, pathout_data_dd, pathout_plots_dd_mean_pair_sync, pathout_plots_dd_sigma_chi, ...
		pathout_plots_dd_pair_sync, pathout_plots_dd_synchrony, pathout_plots_pract_ce_control_vars)
	
%} 

%% scatter plots for sigma met mean & sigma chi mean, with fixed A, and varying beta

% loop over all values of npoints 

%{

get_all_kuramoto_met_chi_scatterplots(network, A_vec, beta_vec, all_npoints, x_label_chi_met, y_label_chi_met, ...
		pathout_data_sync, pathout_plots_sigma_chi, pathout_plots_sigma_met)
	
%}
