function get_ssi_from_var(n_nodes, time_lag_for_model, same_time_coupling_matrix, ...
		full_coupling_matrix, spectral_radius, rmi, noise_corr) 
	
	decay_param		= 1;			% VAR coefficients decay parameter
	frequency_res	= [];			% frequency resolution (empty for automatic)
	moddir		= tempdir;		% model directory
	modname		='sim_model';	% model filename root
	gvdisp		= [];			% GraphViz display? Empty for no action, true to display, false to just generate files)
	gvprog		= 'neato';		% GraphViz program/format (also try 'neato', 'fdp')
	nsics			= 0;			% number of samples for spectral integration check (0 for no check)
	mseed			= 0;			% model random seed (0 to use current rng state)
						
	rstate = rng_seed(mseed);
	
	random_var_model = var_rand(full_coupling_matrix,[],spectral_radius,decay_param);
	
	mdescript = sprintf('%d-variable VAR(%d)',n_nodes,time_lag_for_model);
	causal_graph = var_to_pwcgc(full_coupling_matrix,noise_corr);						% causal graph
	[higher_dim_coupling_matrix,observed_state_var,kalman_var] = var_to_ss(full_coupling_matrix);   % equivalent ISS model
	if isempty(frequency_res)
		[frequency_res,ierr] = var2fres(full_coupling_matrix,noise_corr);
	end
	modcomp = time_lag_for_model*n_nodes*n_nodes + (n_nodes*(n_nodes+1))/2;					% model complexity
	
	rng_restore(rstate);
	
	if nsics > 0
		assert(exist('mdim','var'),'For spectral accuracy check, must supply macro dimension ''mdim''');
		derr = dds_check(A,C,K,H,mdim,nsics);									% spectral integration check
		if derr > 1e-12,
			fprintf(2,'WARNING: spectral DD calculation may be inaccurate!\n\n');
		end
	end
	
	% model info
	
	fprintf('\n---------------------------------------\n');
	fprintf('Model            : %s\n',mdescript);
	fprintf('---------------------------------------\n');
	fprintf('Dimension        : %d\n',n_nodes);
	fprintf('Complexity       : %d\n',modcomp);
	fprintf('Freq. resolution : %d\n',frequency_res);
	fprintf('---------------------------------------\n\n');
	
	if ~isempty(gvdisp)
		eweight = causal_graph/nanmax(causal_graph(:));
		gfile = fullfile(tempdir,'sim_model');
		wgraph2dot(n_nodes,eweight,gfile,[],gvprog,gvdisp);
		fprintf('\n');
	end
	
	modfile = [fullfile(moddir,modname) '.mat'];
	fprintf('*** saving model in ''%s''... ',modfile);
	
	save(modfile,'mdescript','same_time_coupling_matrix','n_nodes','time_lag_for_model','spectral_radius','rmi','noise_corr','full_coupling_matrix','higher_dim_coupling_matrix','observed_state_var','kalman_var','causal_graph','frequency_res','modcomp');
	fprintf('done\n\n');
end 
