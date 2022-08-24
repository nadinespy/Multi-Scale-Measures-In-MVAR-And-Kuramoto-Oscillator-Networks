function get_km_variables(network, time_series_length, model_params, coupling_matrices, ...
		pathout_data_sim_time_series, pathout_data_sync);
	
	% Reminder: model_params = {A, beta, intra_comm_size, n_communities};	% model parameters for kuramoto oscillators;
													% must be in that order
	model_param1 = model_params{1};
	model_param2 = model_params{2};
	model_param3 = model_params{3};
	model_param4 = model_params{4};
	
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

	for q = 1:length(time_series_length);
		time_series_length_str = num2str(time_series_length(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			model_param1_str = param2str(model_param1(i));
		
			for j = 1:length(model_param2)
				model_param2_str = param2str(model_param2(j));
				
				fprintf('get_all_kuramoto_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d\n', q, i, j);
		
				% simulate km oscillators using Shanahan's code
				[phase, chi, sync] = sim_km_oscillators(time_series_length(q), model_param2(j), model_param3, ...
					model_param4, coupling_matrix);

% 				% simulate km oscillators using Lionel's code
% 				w = ones(length(coupling_matrix),1);	                  % identical natural frequencies
% 				h = 0.05;									% Runge-Kutta method step size
% 				initial_phases = rand(1,length(coupling_matrix))*2*pi-pi;	% initial phases of oscillators
% 				dt = 0.01;									% integration time increment
% 				
% 				% to run that function, see instructions here: https://github.com/lcbarnett/kuramoto/blob/main/README.md,
% 				% or run [make -C C && make -C Matlab] in the command line when in the kuramoto directory
% 				[phase,order_param_mag,order_param_phase,sim_time,int_steps] = kuramoto(length(coupling_matrix),w,coupling_matrix,initial_phases,time_series_length(q),dt,1);
% 														% might get warning that LBFGS solver 
% 														% (Limited Broyden–Fletcher–Goldfarb–Shanno algorithm) 
% 														% failed
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


% 				load([pathout_data_sim_time_series network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
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
				save([pathout_data_sim_time_series network '_rica' num2str(n_features1) '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				sum_rica_phase = zeros(1, time_series_length(q));
				for k = 1:(size(rica_phase,1));
					sum_rica_phase = sum_rica_phase + rica_phase(k,:);
				end 
				S.(['sum_rica' num2str(n_features1) '_phase']) = sum_rica_phase;
				save([pathout_data_sim_time_series network '_sum_rica' num2str(n_features1) '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				n_features2 = size(phase,1)/2;
				reconstruction_ica = rica(phase', n_features2); %,'IterationLimit',100);
				rica_phase = (phase' * reconstruction_ica.TransformWeights)';
				S.(['rica' num2str(n_features2) '_phase']) = rica_phase;
				save([pathout_data_sim_time_series network '_rica' num2str(n_features2) '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'-struct', 'S');
				
				clear S;
				
				sum_rica_phase = zeros(1, time_series_length(q));
				for k = 1:(size(rica_phase,1));
					sum_rica_phase = sum_rica_phase + rica_phase(k,:);
				end 
				S.(['sum_rica' num2str(n_features2) '_phase']) = sum_rica_phase;
				save([pathout_data_sim_time_series network '_sum_rica' num2str(n_features2) '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
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
				save([pathout_data_sim_time_series network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'phase');
				save([pathout_data_sim_time_series network '_chi_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'chi');
				save([pathout_data_sim_time_series network '_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'sync');
				save([pathout_data_sim_time_series network '_p_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'p_sync');
				save([pathout_data_sim_time_series network '_mp_sync_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'mp_sync');
				save([pathout_data_sim_time_series network '_raw_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'raw');
				save([pathout_data_sim_time_series network '_sum_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'sum_phase');
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout_data_sync network '_synchronies_'  model_param1_str '_' time_series_length_str '.mat'], ...
				'synchronies');
 			
			clear synchrony;
			clear synchronies;
			
		end
	end
end 

