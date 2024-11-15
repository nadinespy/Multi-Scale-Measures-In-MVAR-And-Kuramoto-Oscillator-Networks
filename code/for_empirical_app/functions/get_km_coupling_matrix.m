function km_coupling_matrix = get_km_coupling_matrix(intra_comm_size, n_communities, A); 
% get_km_coupling () - generates a coupling matrix of size 
% ([intra_comm_size] x [n_communities]) x ([intra_comm_size] x [n_communities])
% and given coupling ratio for difference between intra- and inter-community
% coupling (specified by [A]).
%
% Example: get_km_coupling_matrix(intra_comm_size, n_communities, A)
%
% INPUTS - required:
%    intra_comm_size        -             number of inter-community nodes (double)
%    n_communities          -             number of communities (double) 
%
% OUTPUT:
%    km_coupling_matrix     -             coupling matrix of size 
%                                         ([intra_comm_size] x [n_communities]) x 
%                                         ([intra_comm_size] x [n_communities]) x [A]
% 
% Reference:
%   Shanahan (2010). Metastable chimera states in community-structured oscillator networks.
%
% Code from https://github.com/mpshanahan/metastable-chimera
% and ddapted by Nadine Spychala, Dec 2022
	
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
				km_coupling_matrix(i, j) = k;
				km_coupling_matrix(j, i) = k;
			end
		end
	end
end

end 
