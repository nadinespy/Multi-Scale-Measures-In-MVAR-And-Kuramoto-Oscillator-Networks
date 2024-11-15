% Function to calculate mutual information assuming Gaussianity
% Function to calculate mutual information assuming Gaussianity with PCA dimensionality reduction
function mi = calculateGaussianMutualInformation(x, y)
    % Perform PCA for dimensionality reduction
    [~, score_x] = pca(x, 'NumComponents', 10); % You can adjust the number of components according to your data
    [~, score_y] = pca(y, 'NumComponents', 10); % You can adjust the number of components according to your data
    
    % Compute the covariance matrix
    covMatrix = cov([score_x, score_y]);
    
    % Extract the submatrices of covMatrix
    cov_xx = covMatrix(1:size(score_x, 2), 1:size(score_x, 2));
    cov_yy = covMatrix(size(score_x, 2)+1:end, size(score_x, 2)+1:end);
    cov_xy = covMatrix(1:size(score_x, 2), size(score_x, 2)+1:end);
    
    % Calculate the determinant of the covariance matrices
    det_cov_xx = det(cov_xx);
    det_cov_yy = det(cov_yy);
    det_cov = det(covMatrix);
    
    % Compute the mutual information
    mi = 0.5 * log(det_cov_xx * det_cov_yy / det_cov);
end