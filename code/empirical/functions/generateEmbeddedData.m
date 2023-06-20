% Function to generate embedded data matrix
function embeddedData = generateEmbeddedData(timeSeries, m, tau)
    [n, d] = size(timeSeries); % n: number of time points, d: number of variables
    numEmbeddedPoints = n - (m - 1) * tau;
    embeddedData = zeros(numEmbeddedPoints, m * d);
    
    for i = 1:numEmbeddedPoints
        idx = (i - 1) + (0:(m - 1)) * tau + 1;
        embeddedData(i, :) = reshape(timeSeries(idx, :), 1, []);
    end
end

