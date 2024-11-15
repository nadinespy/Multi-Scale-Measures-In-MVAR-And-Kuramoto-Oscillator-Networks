function mvar_coupling_matrix = get_mvar_coupling_matrix(n_nodes, coupling); 
% get_mvar_coupling () - generates a coupling matrix for [n_nodes] and given 
% coupling value (specified by [coupling], all nodes are connected to each 
% other with the same coupling)
%
% Example: get_mvar_coupling_matrix(n_nodes, coupling)
%
% INPUTS - required:
%    n_nodes                -             number of nodes (double) 
%    coupling               -             coupling value (double) 
%
% OUTPUT:
%    mvar_coupling_matrix   -             coupling matrix of size 
%                                         [n_nodes] x [n_nodes] and
%                                         [coupling] as values
                                                 
	mvar_coupling_matrix = coupling*ones(n_nodes);
end 