function get_km_variables(network, model_sim_params, measure_params, coupling_matrices, ...
		pathout);
% get_km_variables() - generates micro and macro variables in Kuramoto oscillators.
% 
% Takes as inputs, amongst others, arrays with measure parameter values common to 
% all emergence calculations ([measure_params], required), model parameters used 
% for simulation ([model_sim_params], required), and coupling matrices (coupling_matrices,
% required). It then loops over time-lengths and model parameters.
%
% Example: get_km_variables(network, model_sim_params, measure_params, ...
%          coupling_matrices, pathout)
% 
% INPUTS - required:	
%    network                -             char array
%    model_sim_params	    -		      1x1 struct with fields 
%
%							'A': array with doubles 
%                                         'beta': array with doubles	
%							'intra_comm_size': double
%							'n_communities': double
%
%    measure_params         -             1x1 struct with fields
%
%							'measures': cell array with chars
%							'methods': cell array with chars
%							'time_lags': double array
%							'time_lengths': double array
%							'kraskov_params': double array
%							'disc_methods': cell arrays with chars
%							'bins': int array
%
%    coupling_matrices      -             3D double matrix of size 
%							[intra_comm_size x n_communities]
%							 x [intra_comm_size x n_communities]
%							 x [length(A)]
%
%    pathout                -             1x2 struct with field 
%							indicating path to output 
%							for simulated time-series, 
%							and synchronies
%
% OUTPUTS: 
%    micro and macro        -			saved in pathout
%    variables	 
%							

	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'model_sim_params', @isstruct);
	addRequired(p,'measure_params', @isstruct);
	addRequired(p,'coupling_matrices', @isdouble);
	addRequired(p,'pathout', @isstruct);
	
	parse(p, network, model_sim_params, measure_params, coupling_matrices, pathout);
	
	network				= p.Results.network;
	model_sim_params			= p.Results.model_sim_params;
	measure_params			= p.Results.measure_params;
	coupling_matrices			= p.Results.coupling_matrices;
	pathout				= p.Results.pathout;
	
	% extract cell arrays from structs
	pathout_fieldnames		= fieldnames(pathout);
	model_sim_params_fieldnames	= fieldnames(model_sim_params);
	
	pathout1				= pathout.(pathout_fieldnames{1});
	pathout2				= pathout.(pathout_fieldnames{2});
	model_params1			= model_sim_params.(model_sim_params_fieldnames{1});
	model_params2			= model_sim_params.(model_sim_params_fieldnames{2});
	model_params3			= model_sim_params.(model_sim_params_fieldnames{3});
	model_params4			= model_sim_params.(model_sim_params_fieldnames{4});
	time_lengths			= measure_params.time_lengths;
	
	% sim_kuramoto_oscillators() obtains variables 'phase', 'sigma_chi', and 'synchrony';
	% synchronies for different values of beta and given value of A are stored and saved in 'synchronies';
	% 'grand_mean_pair_sync' and 'raw_values' are derived using 'synchrony' and 'phase', respectively

	% SIMULATE KURAMOTO OSCILLATORS: outputs of sim_kuramoto_oscillators() are
	%	- phases				(MICRO)
	%	- chimera-index			(MACRO)
	%	- synchronies			(MICRO)

	% GET FURTHER MICRO AND MACRO VARIABLES:
	%	- raw signal (cos(phase))	(MICRO)
	%	- average pairwise synchrony	(MACRO)

	for q = 1:length(time_lengths);
		time_length_str = num2str(time_lengths(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			model_param1_str = param2str(model_params1(i));
		
			for j = 1:length(model_params2)
				model_param2_str = param2str(model_params2(j));
				
				fprintf('get_kuramoto_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d\n', q, i, j);
		
				% simulate km oscillators using Shanahan's code
				[phase, chi, sync] = sim_km_oscillators(time_lengths(q), model_params2(j), model_params3, ...
					model_params4, coupling_matrix);

% 				% simulate km oscillators using Lionel's code
% 				wmean = 0;									% identical natural frequencies
% 				wsdev = pi/7;
% 				h = 0.05;									% Runge-Kutta method step size
% 				initial_phases = rand(1,length(coupling_matrix))*2*pi-pi;	% initial phases of oscillators
% 				dt = 0.01;									% integration time increment
% 				
% 				wmean = 0;		% oscillator frequencies mean
% 				wsdev = pi/7;	% oscillator frequencies std. dev.
% 				wseed = [];		% oscillator frequencies random seed (empty for no seeding)
% 				Kmean = 0.8;	% oscillator coupling constants mean
% 				Ksdev = Kmean/6;	% oscillator coupling constants std. dev.
% 				Kseed = [];		% oscillator coupling constants random seed (empty for no seeding)
% 				a     = [];		% oscillator phase lag constant
% 				hseed = [];		% oscillator initial phases random seed (empty for no seeding)
% 				T     = 200;	% simulation time
% 				dt    = 0.01;	% integration time increment
% 				nmean = 0.05;	% oscillator input noise magnitude mean (zero for no noise)
% 				nsdev = nmean/3;	% oscillator input noise magnitude std. dev.
% 				nseed = [];		% oscillator input noise magnitude random seed (empty for no seeding)
% 				Iseed = [];		% oscillator input noise random seed (empty for no seeding)
% 				
% 				kuramoto_demo();
% 				
% 				w = wmean + wsdev*randn(length(coupling_matrix),1); % oscillator frequencies normally distributed with mean wmean and std. dev wsdev
% 				% to run that function, see instructions here: https://github.com/lcbarnett/kuramoto/blob/main/README.md,
% 				% or run [make -C C && make -C Matlab] in the command line when in the kuramoto directory
% 				[phase,order_param_mag,order_param_phase,sim_time,int_steps] = kuramoto(length(coupling_matrix),w,coupling_matrix,initial_phases,time_lengths(q),dt,1);
% 														% might get warning that LBFGS solver 
% 														% (Limited Broyden–Fletcher–Goldfarb–Shanno algorithm) 
% 														% failed
% 														
% 				% K is the coupling matrix!
% 				[h2,r2,psi2] = kuramoto(length(coupling_matrix),w,K,a,n,dt,I,'RK4');
% 				
% 				% phase			phase variable (unwrapped)           (N x n matrix)
% 				% order_param_mag		order parameter magnitude            (row vector of length n)
% 				% order_param_phase	order parameter phase (wrapped)      (row vector of length n)
% 				% sim_time			simulation time (possibly adjusted)  (positive double)
% 				% int_steps			integration time steps               (positive integer)

				
				% store synchronies for a given A, and across beta;
				% rows: betas; columns: communities; 3rd dimension: time-points
				synchronies(j,:,:) = sync;
			
				% MICRO VARIABLES:
				%	- phases,
				%	- raw signal,
				%	- synchronies
				%	- pairwise synchrony
				%	- components of phases (same number of components as phases)
				%	- components of phases (half the number of phases)


% 				load([pathout1 network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
% 					'phase');
				
				% raw signal
				raw = cos(phase);

				% do reconstruction ICA with phase & raw signal (number of features to be extracted as high 
				% as number of micro variables): first get weights for each variable and each feature (output of rica()
				% will be a matrix of size [size(phase,2) * number of features], then multiply this matrix with input matrix 
				% to get time-series of independent components/projection of each data point in the component space 
				n_features1 = size(phase,1);
				reconstruction_ica = rica(phase', n_features1); %,'IterationLimit',100);
				rica_phase = (phase' * reconstruction_ica.TransformWeights)';
				S.(['rica' num2str(n_features1) '_phase']) = rica_phase;
				save([pathout1 network '_rica' num2str(n_features1) '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				sum_rica_phase = zeros(1, time_series_length(q));
				for k = 1:(size(rica_phase,1));
					sum_rica_phase = sum_rica_phase + rica_phase(k,:);
				end 
				S.(['sum_rica' num2str(n_features1) '_phase']) = sum_rica_phase;
				save([pathout1 network '_sum_rica' num2str(n_features1) '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				n_features2 = size(phase,1)/2;
				reconstruction_ica = rica(phase', n_features2); %,'IterationLimit',100);
				rica_phase = (phase' * reconstruction_ica.TransformWeights)';
				S.(['rica' num2str(n_features2) '_phase']) = rica_phase;
				save([pathout1 network '_rica' num2str(n_features2) '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				sum_rica_phase = zeros(1, time_series_length(q));
				for k = 1:(size(rica_phase,1));
					sum_rica_phase = sum_rica_phase + rica_phase(k,:);
				end 
				S.(['sum_rica' num2str(n_features2) '_phase']) = sum_rica_phase;
				save([pathout1 network '_sum_rica' num2str(n_features2) '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				% MACRO VARIABLES for practical measures for causal emergence:
				%	- synchronies
				%	- pairwise synchrony
				%	- variance of synchronies (sigma_chi)
				%	- mean pairwise synchrony between communities (mean_pair_sync)
				%	- sum of phases
				%	- sum of components of phases (same number of components as phases)
				%	- sum of components of phases (half the number of phases) 				
 				
				% global mean pairwise synchrony
				[p_sync, mp_sync] = get_kuramoto_pair_sync(model_param4, sync, time_series_length(q));

				% sum over all phases
				sum_phase = zeros(1, time_series_length(q));
				for k = 1:(size(phase,1));
					sum_phase = sum_phase + phase(k,:);
				end 

				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				save([pathout1 network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'phase');
				save([pathout1 network '_chi_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'chi');
				save([pathout1 network '_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'sync');
				save([pathout1 network '_p_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'p_sync');
				save([pathout1 network '_mp_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'mp_sync');
				save([pathout1 network '_raw_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'raw');
				save([pathout1 network '_sum_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'sum_phase');
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout2 network '_synchronies_'  model_param1_str '_' time_length_str '.mat'], ...
				'synchronies');
 			
			clear synchrony;
			clear synchronies;
			
		end
	end
end 

