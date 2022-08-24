function coupling_matrices = get_km_coupling_matrices(get_km_coupling_matrix, model_params)

	% get coupling matrices
	for i = 1:length(model_params{1})
		coupling_matrices(:,:,i) = get_km_coupling_matrix(model_params{2}, model_params{3}, model_params{1}(i));
	end

end
