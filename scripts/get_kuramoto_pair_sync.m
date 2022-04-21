function [pair_sync, mean_pair_sync] = get_kuramoto_pair_sync(n_communities, synchrony, npoints)
	
	% this function calu
	pair_sync = [];
	for h = 1:n_communities;
		for g = 1:n_communities;
			pair_sync = [pair_sync; abs((synchrony(h,:)+synchrony(g,:))/2)];
		end
	end
	
	mean_pair_sync = mean(pair_sync,1);
end 
