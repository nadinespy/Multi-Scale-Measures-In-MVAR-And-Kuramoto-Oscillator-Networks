function [phase,full_system_sync,order_param_phase] = sim_km_oscillators_with_prep(coupling_matrix, phase_lag, time_length);
				
	% simulate km oscillators using Lionel's code
	wmean  = 0;			% identical natural frequencies
	wsdev  = pi/7;		% oscillator frequencies std. dev.
	h	 = 0.05;		% Runge-Kutta method step size
	dt	 = 0.01;		% integration time increment
	wmean  = 0;			% oscillator frequencies mean
	wsdev  = pi/7;		% oscillator frequencies std. dev.
	hseed  = 1;			% oscillator initial phases random seed (empty for no seeding)
	dt	 = 0.01;		% integration time increment
	nmean	 = 0.05;		% oscillator input noise magnitude mean (zero for no noise)
	nsdev	 = nmean/3;		% oscillator input noise magnitude std. dev.
	nseed	 = 1;			% oscillator input noise magnitude random seed (empty for no seeding)
	Iseed	 = 1;			% oscillator input noise random seed (empty for no seeding)
	hifreq = 0.01;		% hi-pass frequency, or 0 for no high-pass
	pdto   = 10;		% polynomial detrend order
	
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
	
	[phase,full_system_sync,order_param_phase] = kuramoto(length(coupling_matrix), w, coupling_matrix, ...
		model_params2(j),n_points,dt,I, 'Euler');
	
	% plot oscillators before preprocessing for full data length
	figure('Position', [10 10 1000 700]);
	for p = 1:length(coupling_matrix)
		subplot(6,4,p);
		plot(phase(p,:));
		subplot(6,4,p+length(coupling_matrix));
		histogram(phase(p,:));
	end
	
	sgtitle('time-series and histograms of oscillators', 'FontSize', 20);
	
	yes_or_no = input('Look at chunks of unprocessed phase? If yes, type 1, if no, type 0, then enter.');
	close all;
	
	% plot data chunks of oscillators before preprocessing
	if yes_or_no == 1
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
	end
	
	% detrend linearly (& high-pass filter)
	if hifreq > 0
		fprintf('\nHigh-pass filter at %g/pi\n',hifreq);
		ppstr = sprintf('High-pass filtered at %g/\\pi',hifreq);
		[fb,fa] = butter(2,hifreq/pi,'high');
		phase = detrend(phase,1);
		phase = filtfilt(fb,fa,phase);
	else
		fprintf('\nPolynomial detrend at order %d\n',pdto);
		ppstr = sprintf('%dth-order polynomial detrend',pdto);
		phase = detrend(phase,pdto);
	end
	
	% plot time-series and histogram for full data length
	figure('Position', [10 10 1100 800]);
	for p = 1:length(coupling_matrix)
		subplot(8,2,p);
		histogram(phase(p,:));
	end
	subplot(8,2,length(coupling_matrix)+1);
	for p = 1:length(coupling_matrix)
		plot(phase(p,:));
		hold on
	end
	subplot(8,2,length(coupling_matrix)+2);
	plot(full_system_sync);
	title('order parameter magnitude');
	subplot(8,2,length(coupling_matrix)+3);
	histogram(full_system_sync);
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
				plot(phase(p,chunk_length*(k-1)+1:chunk_length*k));
				subplot(6,4,p+length(coupling_matrix));
				histogram(phase(p,chunk_length*(k-1)+1:chunk_length*k));
			end
			
			sgtitle(['time-series and histograms of preprocessed oscillators for data chunk: ' ...
				num2str(chunk_length*(k-1)+1) ':' num2str(chunk_length*k)], 'FontSize', 20);
			
			yes_or_no = input('Keep or reject that chunk of data? (Type 0 for no, 1 for yes, then enter.) ');
			all_yes_or_no = [all_yes_or_no; yes_or_no];
			
			close all;
			
			figure;
			plot(full_system_sync(chunk_length*(k-1)+1:chunk_length*k));
			input('Had a look at full_system_sync? If yes, type enter. ');
			
			close all;
			
		end
		
		% exclude data chunks for which a 'no' has been recorded in [all_yes_or_no]
		new_full_system_sync = [];
		new_phase = [];
		for k = 1:n_chunks;
			if all_yes_or_no(k) == 1;
				phase_chunk = phase(:,chunk_length*(k-1)+1:chunk_length*k);
				full_system_sync_chunk = full_system_sync(chunk_length*(k-1)+1:chunk_length*k);
				new_phase = [phase, phase_chunk];
				new_full_system_sync = [new_full_system_sync, full_system_sync_chunk];
			end
		end
		full_system_sync = new_full_system_sync;
		phase = new_phase;
	end
end