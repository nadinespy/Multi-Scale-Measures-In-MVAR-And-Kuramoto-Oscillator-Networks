% Function to calculate mutual information using k-nearest neighbors
function mi = calculate_knn_mi(x, y)
    numPoints = length(x);
    k = ceil(sqrt(numPoints));
    
    % Find k-nearest neighbors
    [~, dist] = knnsearch([x y], [x y], 'K', k + 1);
    dist = dist(:, 2:end); % Exclude the first column (distance to self)
    
    % Calculate mutual information using k-nearest neighbors
    mi = mean(psi(k) + psi(numPoints) - mean(psi(sum(dist, 2))));
end