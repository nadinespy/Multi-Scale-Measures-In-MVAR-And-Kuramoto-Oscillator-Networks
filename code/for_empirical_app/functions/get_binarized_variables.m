function bin_matrix = get_binarized_variables(matrix, bin_threshold)
	
	for k = 1:size(matrix, 2);
		for l = 1:size(matrix,1);
			if matrix(l,k) >= bin_threshold;
				bin_matrix(l,k) = 1;
			else bin_matrix(l,k) = 0;
			end
		end
	end
			
end 