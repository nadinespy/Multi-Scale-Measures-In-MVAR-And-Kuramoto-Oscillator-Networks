function grand_mean_pair_sync = get_kuramoto_grand_mean_pair_sync(n_communities, synchrony, npoints)
	grand_mean_pair_sync = zeros(1, npoints);
	
	for h = 1:n_communities;
		mean_pair_sync_temp = zeros(1, npoints);
		for g = 1:n_communities;
			mean_pair_sync_temp = mean_pair_sync_temp + abs((synchrony(h,:)+synchrony(g,:))/2);
		end
		mean_pair_sync_temp = mean_pair_sync_temp/n_communities;
		grand_mean_pair_sync = grand_mean_pair_sync + mean_pair_sync_temp;
	end
	grand_mean_pair_sync = grand_mean_pair_sync/n_communities;
end 