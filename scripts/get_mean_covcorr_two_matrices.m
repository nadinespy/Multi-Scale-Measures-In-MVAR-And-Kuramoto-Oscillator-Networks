function [mean_cov, mean_corr] = get_mean_covcorr_two_matrices(matrixA, matrixB); 
	
	store_covariance = [];
	store_correlation = [];
	
	for b = 1:size(matrixA, 2);
		for c = 1:size(matrixB,2);
			covariance = cov(matrixA(:,b), matrixB(:,c));
			store_covariance = [store_covariance;covariance(1,2)];
			
			correlation = corrcov(covariance);
			store_correlation = [store_correlation;correlation(1,2)];
		end
	end
	
	mean_cov = mean(store_covariance);
	mean_corr = mean(store_correlation);
	
end 