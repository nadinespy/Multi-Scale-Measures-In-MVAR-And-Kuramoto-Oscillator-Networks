function coupling_matrices = get_km_coupling_matrices(get_coupling_matrix, model_sim_params)

	% Function description: get_km_coupling_matrices() takes as inputs  
	% a function (required) for generating one coupling matrix specifically 
	% for Kuramoto data, and model parameters for simulating the model 
	% ([model_sim_params], required).
	% It then loops over the values of A in model_sim_params. 
	
	% Inputs:	
	%
	% Required:	get_coupling_matrix		function name of the function 
	%							that generates a coupling matrix
	%
	%		model_sim_params			1x1 struct with fields 
	%							A (vector with doubles), beta, 
	%							(vector with doubles), 
	%							intra_comm_size (int),
	%							n_communities (int)
	%
	% Outputs: coupling_matrices			3D double matrix of size 
	%							[intra_comm_size*n_communities
	%							 x intra_comm_size*n_communities
	%							 x length(A)]
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'get_coupling_matrix', @isfunction);
	addRequired(p,'model_sim_params', @isstruct);

	parse(p, get_coupling_matrix, model_sim_params);
	
	model_sim_params			= p.Results.model_sim_params;
	get_coupling_matrix		= p.Results.get_coupling_matrix;
	
	model_sim_params_fieldnames	= fieldnames(model_sim_params);
	model_params1			= model_sim_params.(model_sim_params_fieldnames{1});
	model_params3			= model_sim_params.(model_sim_params_fieldnames{3});
	model_params4			= model_sim_params.(model_sim_params_fieldnames{4});

	
	% get coupling matrices
	for i = 1:length(model_params1)
		coupling_matrices(:,:,i) = get_coupling_matrix(model_params3, model_params4, model_params1(i));
	end

end
