function shuffled_rows = shuffle_rows(matrix)
	
	for k = 1:size(matrix,1)
		idx = randperm(size(matrix,2));
		shuffled_rows(k,idx) = matrix(k,:) ; %corrected
	end
	
end 