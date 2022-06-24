function [quantiles_matrix, quantilized_matrix] = get_quantilized_variables(matrix, quantile_number)
	
	quantiles_matrix = quantile(matrix,quantile_number,2);
	states = [1:quantile_number];
	
	for k = 1:size(matrix, 1);
		for q = 1:size(matrix, 2);

			for f = 1:length(states);
				
				if matrix(k,q) <= quantiles_matrix(k,f)
					quantilized_matrix(k,q) = states(f);
					break
				end 
				
				if f == length(states);
					quantilized_matrix(k,q) = states(f)+1;
				end 
			end 
				
		end 
	end
			
end