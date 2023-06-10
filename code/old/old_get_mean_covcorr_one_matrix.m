function [mean_cov, mean_corr] = get_mean_covcorr_one_matrix(matrix);
	cov_matrix = cov(matrix);
	mean_cov = mean(nonzeros(tril(cov_matrix,-1)), 'all');
			
	corr_matrix = corrcov(cov_matrix);
	mean_corr = mean(nonzeros(tril(corr_matrix,-1)), 'all');
end 
			