



seed_osc_freq = 0; % oscillator frequencies random seed (zero for no seeding)
seed_coup_ratio_mean = 0; % Shanahan connectivity parameter random seed (zero for no seeding)
seed_noise_corr = 0; % oscillator input noise correlation random seed (zero for no seeding)
seed_noise_gen = 0; % oscillator input noise generation random seed (zero for no seeding)
seed_init_phase = 0; % oscillator initial phases random seed (zero for no seeding)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set up oscillator system %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_rmi = 2.^linspace(rmi_range(1), rmi_range(2), n_rmi);		% rmi values (exponential scale)
all_phase_lag = pi/2 - linspace(phase_lag_range(1), ...		% phase lag
	phase_lag_range(2), n_phase_lag); 
all_coup_ratio_mean = linspace(coup_ratio_mean_range(1), ...	% coupling ratio mean values
	coup_ratio_mean_range(2), n_coup_ratio_mean);

n_nodes = n_comms * n_nodes_per_comm;					% total number of oscillators
osc_freq_dev = osc_freq_rel_dev * osc_freq_mean;			% oscillator frequencies deviation
coup_ratio_mean_dev = coup_ratio_mean_rel_dev * coup_ratio_mean;	% coupling ratio mean deviation

% oscillator natural frequencies

if osc_freq_dev > 0
	seed_state = rng_seed(seed_osc_freq);
	osc_freq = pi*ms_betarnd(osc_freq_mean/pi, ...			% oscillator frequencies Beta 
		osc_freq_dev, n_nodes, 1);					% distributed on [0,pi] with mean 
	rng_restore(seed_state);						% osc_freq_mean and 
else											% std. dev osc_freq_dev
	osc_freq = osc_freq_mean * ones(n_nodes, 1);
end
sampling_freq  = 1/dt_sim_points;						% sampling frequency
nyquist_freq = sampling_freq/2;						% Nyqvist frequency
freq_hz = nyquist_freq * osc_freq/pi;					% frequencies in Hz


% simulation times
dt_equil_points = floor(n_equil_points / dt_sim_points);		% equilibriation time increments
assert(dt_equil_points >= 0, 'Equilibriation time too short,' ...
	' or time increment too large!');
n_equil_points = dt_equil_points * dt_sim_points;			% adjusted equilibriation time

sim_time_increments = floor(n_sim_points/dt_sim_points);		% simulation time increments
assert(sim_time_increments > 0,'simulation time too short,' ...
	' or time increment too large!');
n_sim_points = sim_time_increments * dt_sim_points;			% adjusted simulation time

tot_sim_time_increments = dt_equil_points + sim_time_increments;	% total simulation time increments

%% LOOPS
for h = 1:length(all_coup_ratio_mean)

	% set up community-structured weighted network as in Shanahan (2009)
	seed_state = rng_seed(seed_coup_ratio_mean);
	[weighted_network, communities] = shanahan_network(n_comms, ...
		n_nodes_per_comm, ...
		n_inter_comm_coups, ...
		all_coup_ratio_mean(h), ...
		coup_ratio_mean_dev, ...
		global_coup);
	rng_restore(seed_state);

	for i = 1:length(all_phase_lag) 
		
		for j = 1:length(all_rmi)
		
			% generate noise covariance matrix
			seed_state = rng_seed(seed_noise_corr);
			if adjust_rmi
				noise_cov = sqrt(noise_mag) * corr_rand(n_nodes, -all_rmi(j));
			else
				noise_cov = sqrt(noise_mag) * corr_rand(n_nodes, +all_rmi(j));
			end
			rng_restore(seed_state);

			% generate correlated Gaussian white noise
			seed_state = rng_seed(seed_noise_gen);
			noise_corr = randn(tot_sim_time_increments, n_nodes) * chol(noise_cov); 
			rng_restore(seed_state);

			% initial phases uniform on [-pi,pi]
			seed_state = rng_seed(seed_init_phase);
			initial_phases = pi*(2*rand(1, n_nodes)-1);
			rng_restore(seed_state);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%% Run oscillator system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			% run Kuramoto simulations with specified parameters
			[phases, order_param] = kuramoto(n_nodes, ...
				tot_sim_time_increments, ...
				dt_sim_points, osc_freq, ...
				weighted_network, ...
				all_phase_lag(i), ...
				initial_phases, ...
				noise_corr, ...
				sim_mode);

			% truncate equilibriation
			phases = phases(:, dt_equil_points+1:end);
			order_param = order_param(dt_equil_points+1:end);

			% community order parameters (synchronisation) 
			% % --> macroscopic variable
			order_param_comms = zeros(n_comms, sim_time_increments);
			for k = 1:n_comms
				order_param_comms(k,:) = hypot(mean(cos(phases(communities{k},:))), ...
					mean(sin(phases(communities{k},:))));
			end

			DD = nan(n_local_gc, 1);
			k = 0;
			nfail = 0;
			while k < n_local_gc
				if emp_sample
					phases_sample = phases(:, randi(sim_time_increments));
				else
					phases_sample = pi*(2 * rand(n_nodes,1)-1);
				end
				% locally linearised Kuramoto process
				local_lin_phases = kosgrad(phases_sample, weighted_network, phase_lag);				
				local_spectral_radius = max(real(eig(local_lin_phases)));
				% locally linearised macro variable
				local_lin_order_param_comms = kosmacrograd(phases_sample, communities); 
				% compute Dynamical Dependence (DD) for locally 
				% linearised micro & macro variables
				[DD_temp, err] = vou_to_dd(local_lin_phases, ...
					noise_cov, ...
					local_lin_order_param_comms, ...
					stabilise_care);
				if err == 0 % success!
					k = k+1;
					DD(k) = DD_temp;
					if k == n_local_gc, break; 
					end 
				else
					nfail = nfail+1;
				end
			end

			DDm = mean(DD);
			DDs = std(DD);

		end
	end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Display oscillator system dynamics %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FOR PLOTTING

global have_gvmat
if ~have_gvmat, disnet = false; end

defvar('disnet', false        ); % display connectivity graph?
defvar('gvprog','neato'       ); % GraphViz filter ('dot, 'neato', 'fdp', etc. - see man pages; default: 'dot')

% Spectral parameters

defvar('welch_window_width',  []         ); % Welch nwindow width
defvar('welch_noverlap', []         ); % Welch noverlap
defvar('welch_window_type',  []         ); % Welch window type (e.g., @hann, @flattopwin, @rectwin, etc.; see Matlab docs)
defvar('welch_nfft',     []         ); % Welch nfft = 2*fres

% display connectivity
if disnet
	wgraph2dot(osc_freq/max(osc_freq), weighted_network/max(weighted_network(:)), [tempdir 'kuramoto_sh'], true, gvprog, true);
	return
end

% Detrend phases, normalise by variance

detrended_phases = demean(detrend(phases',1)',true);

% Calculate PSDs
fprintf('*** PSD calculation ... ');

[autospectra, freq] = cpsdx(detrended_phases, sampling_freq, welch_window_width, welch_noverlap, welch_window_type, welch_nfft, true);

fprintf('*** Plotting ... ');

freq_hz = nyquist_freq*osc_freq/pi; % frequencies in Hz

time = linspace(0, n_sim_points, sim_time_increments)';

figure(1); clf;
sgtitle(sprintf('Kuramoto system (Shanahan): n_comms = %d, n_nodes_per_comm = %d, n_inter_comm_coups = %d, \\phase_lag = %g\n', n_comms, n_nodes_per_comm, n_inter_comm_coups, all_phase_lag(i)));

subplot(3,1,1)
semilogx(freq, log(autospectra));
xlim([freq(2), nyquist_freq]);
for i = 1:n_nodes, xline(freq_hz(i)); end
xlabel('Frequency (Hz)');
ylabel('Spectral power (dB)');
title(sprintf('PSD\n'));

subplot(3,1,2)
plot(time, X);
xlabel('Time (seconds)');
ylabel('X', 'Rotation', 0);
title(sprintf('De-trended oscillator phases\n'));

subplot(3,1,3)
plot(time, order_param','k','LineWidth',2);
hold on
plot(time, order_param_comms');
hold off
ylim([0, 1]);
xlabel('Time (seconds)');
ylabel('order parameter','Rotation',0);
title(sprintf('Order parameters\n'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% VOU DD analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

