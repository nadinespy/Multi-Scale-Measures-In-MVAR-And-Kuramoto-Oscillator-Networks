function disc_matrix = get_discretized_variables(matrix, bin_number)
	
	for k = 1:size(matrix, 2);
		[disc_matrix, disc_edges] = discretize(matrix,bin_number);
	end
			
end 