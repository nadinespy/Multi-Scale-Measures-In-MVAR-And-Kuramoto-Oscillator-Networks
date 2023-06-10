function [ allres ] = PlotNetworkSummary()
%%PLOTNETWORKSUMMARY Integration measures in networks with different topologies
%
% This function replicates the results in Figure 5 of Mediano et al. 2019, by
% computing several integrated information measures in several 8-node AR
% networks. Note that, as per sections 2.2.2 and 3 of the paper all measures
% are optimised over even-sized bipartitions of the system.
% 
% Most integration measures are computed with my (Pedro's) private fork of the
% JIDT toolbox by Lizier et al. (https://www.github.com/jlizier/jidt), and the
% user is referred to the JIDT documentation for further information about
% using these estimators. Notably, PhiG is computed with version 5 of Kitazono
% and Oizumi's Phi toolbox
% (https://figshare.com/articles/code/phi_toolbox_zip/3203326/5).
%
% Reference:
%
%   Mediano, P.A.; Seth, A.K.; Barrett, A.B. Measuring Integrated Information:
%   Comparison of Candidate Measures in Theory and Simulation. Entropy 2019,
%   21, 17.
%
% Pedro Mediano, 2018-2020

%% Add paths and load coupling matrices
javaaddpath('infodynamics.jar');
javaaddpath('commons-math3-3.5.jar');
addpath('phi_toolbox');

load('nets.mat');
netNames = fields(nets);
nb_nets = length(netNames);


%% Measures to be calculated and their plot labels
calcNames = {'IntegratedInformation', 'IntegratedInteraction', 'AverageCorrelation', ...
    'DecoderIntegration', 'CausalDensity', 'IntegratedSynergy', 'TimeDelayedMutualInfo',...
    'GeometricIntegration'};

calcLabels = {'\Phi', '\tilde\Phi', '\bar\Sigma', '\Phi^*', '\mathrm{CD}',...
    '\psi', 'I(X_{t-\tau}, X_t)', '\Phi_G'};

nb_calcs = length(calcNames);

% Name template to instantiate JIDT calculators
class_template = 'infodynamics.measures.continuous.gaussian.%sCalculatorGaussian';


%% Loop and compute all measures for all networks
allres = struct();
for calc_idx=1:nb_calcs
    
    calcName = calcNames{calc_idx};

    res = zeros([nb_nets, 1]);
    for i=1:nb_nets

        % Normalise network so that spectral radius is close to but below 1
        net = nets.(netNames{i});
        sr = max(abs(eig(net)));
        net = net./(1.10*sr);
        lS = makeLaggedCovariance(net, 0);
        
        % Compute measure
        if strcmp(calcName, 'AverageCorrelation')
          res(i) = (sum(abs(lS(:))) - sum(abs(diag(lS))))/(16*16-16);
        elseif strcmp(calcName, 'GeometricIntegration')
          res(i) = PhiGMatlab(lS, 'even_bipartitions');
        else
          calc = javaObject(sprintf(class_template, calcName));
          calc.setProperty('PARTITION_SCAN_METHOD', 'EVEN_BIPARTITIONS');
          calc.initialise(size(net, 1));
          calc.setLaggedCovariance(lS);
          res(i) = calc.computeAverageLocalOfObservations();
        end
        
    end

    % Plot and copy results to main struct
    subplot(2,4,calc_idx);
    bar(res);
    title(sprintf('$%s$', calcLabels{calc_idx}), 'Interpreter', 'latex')
    colormap(0.6*ones([1, 3]));
    set(gca, 'XLim', [0.25, 6.75]);
    set(gca, 'XTick', 1:nb_nets, 'XTickLabel', cellstr(char(64+(1:nb_nets))'));
    set(gca, 'FontSize', 20);

    allres.(calcName) = res;
    
end
    
end  % of main function

function [ lS ] = makeLaggedCovariance(A, c)
  %% Auxiliary function to compute analytically the time-lagged covariance
  % matrix of an AR process with coupling matrix A and noise correlation c
  assert(size(A,1) == size(A,2));
  N = size(A, 1);
  eps = c*ones(N);
  eps(1:(N+1):end) = 1;

  S = dlyap(A, eps);
  lS = [S, S*A'; A*S, S];
  lS = 0.5*(lS + lS');
end

