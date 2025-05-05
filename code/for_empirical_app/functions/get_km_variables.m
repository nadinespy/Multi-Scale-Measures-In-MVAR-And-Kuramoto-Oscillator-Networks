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
	
	% sim_kuramoto_oscillators() obtains variables 'phase', 'sigma_chi', and 'community_sync';
	% synchronies for different values of beta and given value of A are stored and saved in 'synchronies';
	% 'grand_mean_pair_sync' and 'raw_values' are derived using 'community_sync' and 'phase', respectively

	% SIMULATE KURAMOTO OSCILLATORS: outputs of sim_kuramoto_oscillators() are
	%	- phases				(MICRO)
	%	- order parameter magnitude   (MACRO)
	%       for full system

	% GET FURTHER MICRO AND MACRO VARIABLES:
	%	- raw signal (cos(phase))	(MICRO)
	%	- synchronies			(MICRO)
	%	- average pairwise synchrony	(MACRO)
	%       of communities
	%	- chimera-index			(MACRO)
	%	- order parameter magnitude   (MACRO)
	%       for communities

	for q = 1:length(time_lengths);
		time_length_str = num2str(time_lengths(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			model_param1_str = param2str(model_params1(i));
		
			for j = 5:length(model_params2)
				model_param2_str = param2str(model_params2(j));
				
				fprintf('get_kuramoto_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d\n', q, i, j);
		
				% -----------------------------------------------------------
				% SIMULATE PHASES & ODER PARAMETER MAGNITUDES
				% -----------------------------------------------------------
				
 				% simulate km oscillators using Shanahan's code
 				% [phase, chi, sync] = sim_km_oscillators(time_lengths(q), model_params2(j), model_params3, ...
 				%	model_params4, coupling_matrix);

				
				
				% to run that function, see instructions here: https://github.com/lcbarnett/kuramoto/blob/main/README.md,
				% or run [make -C C && make -C Matlab] in the command line when in the kuramoto directory;
				% might get warning that LBFGS solver (Limited Broyden–Fletcher–Goldfarb–Shanno algorithm) failed
				
				% phase			phase variable (unwrapped)           (N x n matrix)
				% full_system_sync	order parameter magnitude            (row vector of length n)
				% order_param_phase	order parameter phase (wrapped)      (row vector of length n)
				
				yes_or_no = 0;
				while yes_or_no == 0
					[phase, full_system_sync, order_param_phase] = sim_km_oscillators_with_gc_test(coupling_matrix, model_params2(j), time_lengths(q));
					yes_or_no = input('Keep this KM simulation? If yes, type 1, if no, type 0, then enter.');
					close all;
				end

				% [phase,full_system_sync,order_param_phase] = sim_km_oscillators_with_prep(coupling_matrix, model_params2(j), time_lengths(q));
				
				% -----------------------------------------------------------
				% SYNCHRONIES/ORDER PARAMETER OF COMMUNITIES & CHIMERA-INDEX
				% -----------------------------------------------------------
				community_sync = zeros(length(phase), model_params4);
				chi = zeros(1, length(phase));
				phase = phase';
					
				for c = 1:model_params4
					for e = 1:model_params3
						oscillator = phase(:,(c-1)*model_params3+e);
						community_sync(:,c) = community_sync(:,c)+exp(oscillator*sqrt(-1));	% add all synchrony values of each oscillator belonging to the same community
					end
				end
				
				community_sync = abs(community_sync/model_params3); % take the average
				chi(1,:) = var(community_sync, 0, 2);
							
				community_sync = community_sync';
				phase	    = phase';
				chi	    = chi';
				
% 				% different way of calculating the synchrony/order parameter of 
% 				% communities
% 				% -----------------------------------------------------------
% 				% SYNCHRONIES/ORDER PARAMETER OF COMMUNITIES
% 				% -----------------------------------------------------------
% 				for s = 1:model_params4
% 					temp_community = phase(model_params3*(s-1)+1:model_params3*s,:);
% 
% 					add_exp_community = 0;
% 					for g = 1:(length(coupling_matrix)/model_params4)
% 						exp_community = exp(temp_community(g,:).*sqrt(-1));
% 						add_exp_community = add_exp_community + exp_community;
% 					end
% 					
% 					order_param_mag_community(s,:) = abs((1/(length(coupling_matrix)/model_params4)) * add_exp_community);
% 				end
				
				% store synchronies for a given A, and across beta;
				% rows: betas; columns: communities; 3rd dimension: time-points
				synchronies(j,:,:) = community_sync;
			
				% -----------------------------------------------------------
				% MICRO VARIABLES:
				% -----------------------------------------------------------
				%	- phases,
				%	- raw signal,
				%	- synchronies/order parameter magnitude of communities
				%	- pairwise synchrony of communities
				%	- components of phases (same number of components as phases)
				%	- components of phases (half the number of phases)

 				% load([pathout1 network '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
 				%	'phase');
				
				% -----------------------------------------------------------
				% RAW SIGNAL
				% -----------------------------------------------------------
				raw = cos(phase);
				
				% -----------------------------------------------------------
				% COMPONENTS OF PHASES
				% -----------------------------------------------------------
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
				
				sum_rica_phase = zeros(1, length(rica_phase));
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
				
				sum_rica_phase = zeros(1, length(rica_phase));
				for k = 1:(size(rica_phase,1));
					sum_rica_phase = sum_rica_phase + rica_phase(k,:);
				end 
				S.(['sum_rica' num2str(n_features2) '_phase']) = sum_rica_phase;
				save([pathout1 network '_sum_rica' num2str(n_features2) '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				% -----------------------------------------------------------
				% MACRO VARIABLES:
				% -----------------------------------------------------------
				%	- pairwise synchrony of communities
				%	- variance of synchronies/order parameter magnitudes (sigma_chi)
				%	- mean pairwise synchrony between communities (mean_pair_sync)
				%	- sum of phases
				%	- sum of components of phases (same number of components as phases)
				%	- sum of components of phases (half the number of phases)
				%     - synchrony/order parameter magnitude of full system
				%     - synchronies/order parameter magnitude of communities	
				%     - sum of components of phases (same number of components as phases)
				%	- sum of components of phases (half the number of phases)

 				
				% -----------------------------------------------------------
				% GLOBAL MEAN PAIRWISE SYNCHRONY
				% -----------------------------------------------------------
				[p_sync, mp_sync] = get_kuramoto_pair_sync(model_params4, community_sync, time_lengths(q));

				% -----------------------------------------------------------
				% SUM OVER ALL PHASES
				% -----------------------------------------------------------
				sum_phase = zeros(1, length(phase));
				for k = 1:(size(phase,1));
					sum_phase = sum_phase + phase(k,:);
				end 
				
				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				save([pathout1 network '_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'phase');
				save([pathout1 network '_full_system_sync_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'full_system_sync');
				save([pathout1 network '_chi_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'chi');
				save([pathout1 network '_community_sync_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'community_sync');
				save([pathout1 network '_p_sync_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'p_sync');
				save([pathout1 network '_mp_sync_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'mp_sync');
				save([pathout1 network '_raw_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'raw');
				save([pathout1 network '_sum_phase_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'sum_phase');
				
				j_disp = ['beta index = ',num2str(j)];
				disp(j_disp);
				
				i_disp = ['A index = ',num2str(i)];
				disp(i_disp);
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout2 network '_synchronies_'  model_param1_str '_' time_length_str '.mat'], ...
				'synchronies');
 			
			clear community_sync;
			clear synchronies;
			
		end
	end
end 

