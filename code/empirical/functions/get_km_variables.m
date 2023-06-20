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
	%	- order parameter magnitude   (MACRO)
	%       for full system

	% GET FURTHER MICRO AND MACRO VARIABLES:
	%	- raw signal (cos(phase))	(MICRO)
	%	- synchronies			(MICRO)
	%	- average pairwise synchrony	(MACRO)
	%	- chimera-index			(MACRO)
	%	- order parameter magnitude   (MACRO)
	%       for communities

	for q = 1:length(time_lengths);
		time_length_str = num2str(time_lengths(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			model_param1_str = param2str(model_params1(i));
		
			for j = 1:length(model_params2)
				model_param2_str = param2str(model_params2(j));
				
				fprintf('get_kuramoto_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d\n', q, i, j);
		
				% -----------------------------------------------------------
				% SIMULATE PHASES & ODER PARAMETER MAGNITUDES
				% -----------------------------------------------------------
				
 				% simulate km oscillators using Shanahan's code
 				% [phase, chi, sync] = sim_km_oscillators(time_lengths(q), model_params2(j), model_params3, ...
 				%	model_params4, coupling_matrix);

				% simulate km oscillators using Lionel's code
				wmean  = 0;		% identical natural frequencies
				wsdev  = pi/7;	% oscillator frequencies std. dev.
				h	 = 0.05;	% Runge-Kutta method step size
				dt	 = 0.01;	% integration time increment
				wmean  = 0;		% oscillator frequencies mean
				wsdev  = pi/7;	% oscillator frequencies std. dev.
				hseed  = 1;		% oscillator initial phases random seed (empty for no seeding)
				dt	 = 0.01;	% integration time increment
				nmean	 = 0.05;	% oscillator input noise magnitude mean (zero for no noise)
				nsdev	 = nmean/3;	% oscillator input noise magnitude std. dev.
				nseed	 = 1;		% oscillator input noise magnitude random seed (empty for no seeding)
				Iseed	 = 1;		% oscillator input noise random seed (empty for no seeding)
				hifreq = 0.01;    % hi-pass frequency, or 0 for no high-pass
				pdto   = 10;      % polynomial detrend order
				
				n_points = round(time_lengths(q)/dt);
				assert(n_points > 0,'Simulation time too short, or time increment too large!');
				T = n_points*dt;  % adjusted simulation time
				
				%kuramoto_demo();
				
				if nmean > 0 % with input noise
					if ~isempty(nseed), rstate = rng(nseed); end
					lnv = log(1+nsdev^2/nmean^2);
					% per-oscillator noise magnitudes drawn from log-normal distribution
					nmag = lognrnd(log(nmean)-lnv/2,sqrt(lnv),length(coupling_matrix),1); 
					if ~isempty(nseed), rng(rstate); end
				end
				
				if nmean > 0 % with input noise
					if ~isempty(Iseed), rstate = rng(Iseed); end
					% uncorrelated Gaussian white noise
					I = nmag.*randn(length(coupling_matrix),n_points); 
					if ~isempty(Iseed), rng(rstate); end
				else
					I = zeros(length(coupling_matrix),n_points); % no input
				end
				
				if ~isempty(hseed), rstate = rng(hseed); end
					% initial phases uniform on [-pi,pi]
					I(:,1) = pi*(2*rand(length(coupling_matrix),1)-1);      
				if ~isempty(hseed), rng(rstate); end

				% oscillator frequencies normally distributed with mean wmean and std. dev wsdev
				w = wmean + wsdev*randn(length(coupling_matrix),1); 
				
				% to run that function, see instructions here: https://github.com/lcbarnett/kuramoto/blob/main/README.md,
				% or run [make -C C && make -C Matlab] in the command line when in the kuramoto directory;
				% might get warning that LBFGS solver (Limited Broyden–Fletcher–Goldfarb–Shanno algorithm) failed
				
				% phase			phase variable (unwrapped)           (N x n matrix)
				% order_param_mag		order parameter magnitude            (row vector of length n)
				% order_param_phase	order parameter phase (wrapped)      (row vector of length n)
				[phase,order_param_mag,order_param_phase] = kuramoto(length(coupling_matrix),w,coupling_matrix,model_params2(j),n_points,dt,I, 'Euler');
				
				% plot data chunks of oscillators before preprocessing
				n_chunks = 10;
				chunk_length = n_points/n_chunks;
				for k = 1:n_chunks;
					
					figure('Position', [10 10 1000 700]);
					for p = 1:length(coupling_matrix)
						subplot(6,4,p);
						plot(phase(p,chunk_length*(k-1)+1:chunk_length*k));
						subplot(6,4,p+length(coupling_matrix));
						histogram(phase(p,chunk_length*(k-1)+1:chunk_length*k));
					end
					
					sgtitle(['time-series and histograms of oscillators for data chunk: ' ...
						num2str(chunk_length*(k-1)+1) ':' num2str(chunk_length*k)], 'FontSize', 20);
					
					input('Had a look at this chunk of phase? If yes, type enter.');
					close all;
					
				end
				
				% detrend linearly
				detrend_phase = detrend(phase,1);	
				
				if hifreq > 0
					fprintf('\nHigh-pass filter at %g/pi\n',hifreq);
					ppstr = sprintf('High-pass filtered at %g/\\pi',hifreq);
					[fb,fa] = butter(2,hifreq/pi,'high');
					detrend_phase = filtfilt(fb,fa,detrend_phase);
				else
					fprintf('\nPolynomial detrend at order %d\n',pdto);
					ppstr = sprintf('%dth-order polynomial detrend',pdto);
					detrend_phase = detrend(phase,pdto);
				end
				
				% plot time-series and histogram for full data length
				figure('Position', [10 10 1100 800]);
				for p = 1:length(coupling_matrix) 
					subplot(8,2,p);
					histogram(detrend_phase(p,:));
				end
				subplot(8,2,length(coupling_matrix)+1);
				for p = 1:length(coupling_matrix)
					plot(detrend_phase(p,:));
					hold on
				end
				subplot(8,2,length(coupling_matrix)+2);
				plot(order_param_mag);
				title('order parameter magnitude');
				subplot(8,2,length(coupling_matrix)+3);
				histogram(order_param_mag);
				title('order parameter magnitude');

				sgtitle('time-series and histogram of preprocessed oscillators for full data length', ...
					'FontSize', 20);
				
				look_at_chunks = input('Look at chunks of data? (Type 0 for no, 1 for yes, then enter.) ');
				
				close all;
				
				% visually inspect different chunks of the data for all oscillators after
				% preprocessing, if full-length data hints to non-stationarity/non-Gaussianity; 
				% if given chunk doesn't look good for one or more oscillators, 
				% exclude it from the time-series of all oscillators
				
				% visual inspection; storing rejections/acceptation in [all_yes_or_no]
				if look_at_chunks == 1
					
					n_chunks = 10;
					chunk_length = n_points/n_chunks;
					all_yes_or_no = [];
					for k = 1:n_chunks;
						
						figure('Position', [10 10 1000 700]);
						for p = 1:length(coupling_matrix)
							subplot(6,4,p);
							plot(detrend_phase(p,chunk_length*(k-1)+1:chunk_length*k));
							subplot(6,4,p+length(coupling_matrix));
							histogram(detrend_phase(p,chunk_length*(k-1)+1:chunk_length*k));
						end
						
						sgtitle(['time-series and histograms of preprocessed oscillators for data chunk: ' ...
							num2str(chunk_length*(k-1)+1) ':' num2str(chunk_length*k)], 'FontSize', 20);
						
						yes_or_no = input('Keep or reject that chunk of data? (Type 0 for no, 1 for yes, then enter.) ');
						all_yes_or_no = [all_yes_or_no; yes_or_no];
						
						close all;
						
						figure;
						plot(order_param_mag(chunk_length*(k-1)+1:chunk_length*k));
						input('Had a look at order_param_mag? If yes, type enter. ');
						
						close all;
						
					end
					
					% exclude data chunks for which a 'no' has been recorded in [all_yes_or_no]
					phase = [];
					new_order_param_mag = [];
					for k = 1:n_chunks;
						if all_yes_or_no(k) == 1;
							data_chunk = detrend_phase(:,chunk_length*(k-1)+1:chunk_length*k);
							order_param_mag_chunk = order_param_mag(chunk_length*(k-1)+1:chunk_length*k);
							phase = [phase, data_chunk];
							new_order_param_mag = [new_order_param_mag, order_param_mag_chunk];
						end
					end
					order_param_mag = new_order_param_mag;
					
				end 
				
				% -----------------------------------------------------------
				% SYNCHRONIES/ORDER PARAMETER OF COMMUNITIES & CHIMERA-INDEX
				% -----------------------------------------------------------
				synchrony = zeros(length(phase), model_params4);
				chi = zeros(1, length(phase));
				phase = phase';
				for t = 2:length(phase);
					
					for c = 1:model_params4
						for j = 1:model_params3
							x1 = phase(t,(c-1)*model_params3+j);
							synchrony(t,c) = synchrony(t,c)+exp(x1*sqrt(-1));	% add all synchrony values of each oscillator belonging to the same community
						end
					end
					
					synchrony(t,:) = abs(synchrony(t,:)/model_params3); % take the average
					chi(1,t) = var(synchrony(t,:));
				end 
							
				synchrony = synchrony';
				phase	    = phase';
				chi	    = chi';
				
% 				% different way of calculating the synchrony/order parameter of 
% 				% communities - didn't check yet whether both methods give same results
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
				synchronies(j,:,:) = synchrony;
			
				% -----------------------------------------------------------
				% MICRO VARIABLES:
				% -----------------------------------------------------------
				%	- phases,
				%	- raw signal,
				%	- synchronies/order parameter magnitude of communities
				%	- pairwise synchrony
				%	- components of phases (same number of components as phases)
				%	- components of phases (half the number of phases)

 				% load([pathout1 network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
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
				
				% -----------------------------------------------------------
				% MACRO VARIABLES:
				% -----------------------------------------------------------
				%	- pairwise synchrony
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
				[p_sync, mp_sync] = get_kuramoto_pair_sync(model_param4, sync, time_lengths(q));

				% -----------------------------------------------------------
				% SUM OVER ALL PHASES
				% -----------------------------------------------------------
				sum_phase = zeros(1, time_series_length(q));
				for k = 1:(size(phase,1));
					sum_phase = sum_phase + phase(k,:);
				end 
				
				% variable names consist of  network name + variable name + value of A + value of beta + number of datapoints
				save([pathout1 network '_phase_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'phase');
				save([pathout1 network '_order_param_mag_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'order_param_mag');
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
				save([pathout1 network '_order_param_mag_community_' model_param1_str '_' model_param2_str '_' time_series_length_str '.mat'], ...
					'order_param_mag_community');
			end
		
			% save synchronies; saved filename consists of network name + variable name + value of A + number of datapoints
			save([pathout2 network '_synchronies_'  model_param1_str '_' time_length_str '.mat'], ...
				'synchronies');
 			
			clear synchrony;
			clear synchronies;
			
		end
	end
end 

