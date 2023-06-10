%% phiid-based measures (don't work for this system size)

% some covariance matrices are not positive definite which is why some information-theoretic functions in PhiID 
% can't be computed (in this case, all terms involving t1 and t2, i.e., the indices of first & second target partition):
% h_p1t1t2 = h([p1 t1 t2]);									 
% h_p2t1t2 = h([p2 t1 t2]);
% h_t1t2 = h([t1 t2]);											
% h_p1p2t1t2 = h([p1 p2 t1 t2]);
% where h = @(idx) -log(mvnpdf(sX(idx,:)', mu(idx), S(idx,idx))); 
% h() gives the multivariate entropy, takes as an input the indices of the variables to consider in sX in PhiIDFull()

% one way to check whether a given matrix is positive definite is to see whether p is positive or not - it will be
% zero, if the matrix is positive definite, and positive, if it's not:
% [~,p] = chol(some_matrix)

% calculating information atoms

rng(1);
for i = 1:size(coupling_matrices, 3);
	
	coupling_matrix = coupling_matrices(:,:,i);
	disp(i)
	
	for j = 1:length(beta_vec)
	
		beta = beta_vec(j);
		
		[thetas, sigma_chi, synchrony] = sim_method(coupling_matrix, npoints, beta, intra_comm_size, n_communities);	
		
		% partition indices (arbitrary assigment)
		part1 = linspace(1, 128, 128);				    % indices of partition 1, e.g., in an 8-element system, [1,2,7,8]
		part2 = linspace(129, 256, 128);				  % indices of partition 2, e.g., in an 8-element system, [3,4,5,6]

		% stack data and call full PhiID function
		X1 = thetas(part1,1:end-tau);
		X2 = thetas(part2,1:end-tau);					    % Here, we define X1, X2, Y1, and Y2, so the partitions at time t (X1, X2), and the partitions at time t+1 (Y1, Y2)
            Y1 = thetas(part1,1+tau:end);
		Y2 = thetas(part2,1+tau:end);
		
		X = [X1; X2; Y1; Y2];								% stack all variables from partitions row-wise
		sX = X./repmat(std(X')', [1, npoints - tau]);
		X1 = sX(1:128,:);
		X2 = sX(129:256,:);
		Y1 = sX(257:384,:);
		Y1 = sX(385:512,:);

		% PhiID
		try
			phiid_all_beta_coup_mmi(:,i,j) = struct2array(PhiIDFull(X1, X2, Y1, Y2, 'MMI'))';
		catch 
			phiid_all_beta_coup_mmi(:,i,j) = NaN;
		end
			
		try
			phiid_all_beta_coup_ccs(:,i,j) = struct2array(PhiIDFull(X1, X2, Y1, Y2, 'ccs'))';
		catch 
			phiid_all_beta_coup_ccs(:,i,j) = NaN;
		end
			
	end
	
end 


			% 		% adding error to micro variables
			% 		cov_err  = eye(N, N);
			% 		mu = zeros(N,1);
			% 		E = mvnrnd(mu, cov_err, npoints)';
			% 		raw_signal = raw_signal + E;
			% 		phase = phase + E;
			%
			% 		cov_err  = eye(size(synchrony,1), size(synchrony,1));
			% 		mu = zeros(size(synchrony,1),1);
			% 		E = mvnrnd(mu, cov_err, npoints)';
			% 		synchrony = synchrony + E;