%% integration across RMI and connectivities (where only one coupling changes)

clear all;
close all;
clc;

mvgc_path  = getenv('MVGC2_PATH');
run(fullfile(mvgc_path,'startup'));

cd '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/common'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/analytical/functions'
addpath '/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/code/analytical/scripts/integration_change_one_coupling'
pathout_plots = {'/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/results/plots/2node_mvar/analytical/integration_change_one_coupling/'};
pathout_data_measures = {'/media/nadinespy/NewVolume1/work/current_projects/emergence_complexity_experiments/emergence_complexity_simulations/emergence_complexity_simulations/results/analyses/2node_mvar/analytical/integration_change_one_coupling/'};

n				= 2;		% number of variables in network
time_lag			= 1;		% number of time-lags
n_samples_noise_corrs	= 2;		% number of random samples for noise correlation matrices
seed				= 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed) % seed random number generator

%% create coupling matrices

% variation of residuals' mutual information and matrix norms
all_rmi	= linspace(0.0,1,100);		% array of rmi values 

% two-node matrices of the form [a c; 0 b] for a = 0.8; b = 0.7;
% and c = y where y is one value from [-2:1000:2]
a = 0.8;
b = 0.7;
d = 0;
cr = 2;
c = linspace(-cr,cr,100)';
for i = 1:length(c);

	for i = 1:length(c);
		coupling_matrix = [a c(i); d, b];
		coupling_matrices(:,:,i) = coupling_matrix;
	end

end  
%}

%% loop over coupling matrices and RMI values (including random sampling of noise correlation matrices)

for i = 1:length(coupling_matrices);
	
	for j = 1:length(all_rmi);

			for p = 1:n_samples_noise_corrs;
				
				% define custom error message
				errID = 'DLYAP:NotDefinedOrNotUnique';
				msgtext = 'No or no unique solution to Lyapunov equation.';
				ME = MException(errID,msgtext);
				
				% try to generate noise corrs/coupling matrices that give
				% unique DLYAP solutions, as long as the previous one was neither 
				% unique nor existent
				count = 0;
				err_count = 0;
				while count == err_count 
					
					try
						% append to account for time-lag
						coupling_matrix = [];
						for g = 1:time_lag+1
							coupling_matrix = cat(3,coupling_matrix, coupling_matrices(:,:,i));
						end
					
						% normalize by spectral radius 
						coupling_matrix = specnorm(coupling_matrix,1);
						
						noise_corr  = corr_rand(n,all_rmi(j));
						C		= var_to_autocov(coupling_matrix,noise_corr,time_lag);	% process covariance matrix (solve DLYAP)
						Iv(p)		= sum(log(diag(C(:,:,1))));						% process total variance
						Ig(p)		= log(det(C(:,:,1)));							% process generalised variance (includes covariances)
						I(p)		= Iv(p) - Ig(p);							% process multi-information (nats)
					
					catch % ME.message
						err_count = err_count + 1;
					end
					count = count + 1;
					
					if count == err_count && err_count >= 10
						Iv(p)	= NaN;						
						Ig(p)	= NaN;						
						I(p)	= NaN;
						break
					end
					
				end
					
			end
		
		integration{i, j} = real(I);
		
	end
		
	disp(i)
end

%% save into table

all_c_str = {};
for t = 1:length(c)
	all_c_str{t} = num2str(c(t));
end

all_rmi_str = {};
for e = 1:length(all_rmi)
	all_rmi_str{e} = num2str(all_rmi(e));
end

integration = array2table(integration, 'RowNames', all_c_str, ...
	'VariableNames', all_rmi_str);

save([char(pathout_data_measures), 'integration.mat'],'integration');