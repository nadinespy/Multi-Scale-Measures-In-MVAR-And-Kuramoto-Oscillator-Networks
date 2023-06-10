function [phi_IG, D_A_vec] = phi_G_grad_Ap_copy( x, Cov_E_p, Cov_X, Cov_E, A, partition_cell )


N = size(Cov_X,1);

A_p = zeros(size(A));
idx_st = 0;
for i = 1:length(partition_cell)
    M = partition_cell{i};
    nnz_cell_i = length(M);
    idx_end = nnz_cell_i^2;
    A_p(M,M) = reshape(x(idx_st + (1:idx_end)), [nnz_cell_i, nnz_cell_i]);
    idx_st = idx_st + idx_end;
end

R = [Cov_X Cov_X*A'; A*Cov_X Cov_E+A*Cov_X*A'];
Rd = [Cov_X Cov_X*A_p'; A_p*Cov_X Cov_E_p+A_p*Cov_X*A_p'];
Rd_inv = [inv(Cov_X)+A_p'/Cov_E_p*A_p -A_p'/Cov_E_p; -Cov_E_p\A_p inv(Cov_E_p)];

TR = trace(R*Rd_inv);
phi_IG = 1/2*(-logdet(R) + TR + logdet(Rd) - 2*N);


A_diff = A_p - A;

D_A = 2*Cov_E_p\A_diff*Cov_X;
D_A_vec = zeros(length(x),1);
idx_st = 0;
for i = 1:length(partition_cell)
    M = partition_cell{i};
    nnz_cell_i = length(M);
    idx_end = nnz_cell_i^2;
    D_A_vec(idx_st + (1:idx_end)) = reshape(D_A(M,M), [idx_end, 1]);
    idx_st = idx_st + idx_end;
end


end

