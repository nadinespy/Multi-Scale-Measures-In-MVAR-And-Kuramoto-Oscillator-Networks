function kuramoto_coupling_matrix = get_km_coupling_matrix(intra_comm_size, n_communities, A); 

d0 = intra_comm_size; 	
d1 = intra_comm_size;							% numbers of connections at different community levels
		
N = intra_comm_size*n_communities;					% total number of oscillators: 8
M = n_communities;							% number of lowest level communities (what's that?): 2
	
k1 = (1-A)/2;
k0 = 1-k1;

% build coupling matrix
for i = 1:N
	x1 = mod(ceil(i/intra_comm_size)-1,n_communities)+1;				% community number
	for j = i:N
		if i ~= j										% ignore diagonals
			y1 = mod(ceil(j/intra_comm_size)-1,n_communities)+1;		% community number
			if x1 == y1									% same community
				p = d0/intra_comm_size;
				k = k0;
			else										% different communities
				p = d1/(intra_comm_size*n_communities);
				k = k1;
			end
			if rand < p
				kuramoto_coupling_matrix(i, j) = k;
				kuramoto_coupling_matrix(j, i) = k;
			end
		end
	end
end

end 
