%% GET RIGHT TIME LAG FOR KM OSCILLATORS

tau = 1;	% Initial guess for the time-lag
m = 3;	% Embedding dimension
tauMax = 10; % Maximum time-lag to consider
% Calculate mutual information for each time lag
miValues = zeros(tauMax, 1);

for lag = 1:tauMax
    % Generate the embedded data matrix
    embeddedData = generateEmbeddedData(new_phase, m, lag * tau);
    
    % Extract the variables from the embedded data
    x = embeddedData(:, 1); % First variable
    y = embeddedData(:, 2); % Second variable
    
    % Estimate mutual information using the KDE method
    miValues(lag) = calculate_knn_mi(x, y);
end

% Find the optimal time lag with maximum mutual information
[max_mi, opt_time_lag] = max(miValues);

% Display the mutual information values for different time lags
lagValues = 1:tauMax;
plot(lagValues, miValues, 'o-');
xlabel('Time Lag');
ylabel('Mutual Information');
title('Mutual Information vs. Time Lag');
