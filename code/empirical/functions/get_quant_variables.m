function [quantiles_matrix, quantilized_matrix] = get_quant_variables(matrix, quantile_number)

	states = [1:quantile_number];
	
	% get quantiles in the interval [0,1] based on desired number of quantiles 
	% (e. g., if quantile_number = 1, then all_quantiles = 0.5, or 
	% if quantile_number = 2, then all_quantiles = [0.33, 0.66])
	
	first_quantile = 1/(quantile_number+1);
	rolling_value = 0;
	for j = 1:quantile_number
		rolling_value = rolling_value + first_quantile;
		all_quantiles(j) = rolling_value;
	end 
	
	% get quantiles in the data
	quantiles_matrix = quantile(matrix,all_quantiles,2);
	
	% sort data into states according to quantiles
	for k = 1:size(matrix, 1);
		for q = 1:size(matrix, 2);

			for f = 1:length(states);
				
				if matrix(k,q) <= quantiles_matrix(k,f)
					quantilized_matrix(k,q) = states(f)-1;
					break
				end 
				
				if f == length(states);
					quantilized_matrix(k,q) = states(f);
				end 
			end 
				
		end 
	end
			
end
