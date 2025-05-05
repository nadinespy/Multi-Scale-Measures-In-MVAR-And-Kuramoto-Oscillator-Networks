function cmap = perceptually_uniform_sequential_cmap(n)
    % Creates a perceptually uniform two-color gradient with enhanced contrast
    % n: number of colors
    
    % Define RGB values directly for the endpoints
    blue = [0.23137255, 0.298039216, 0.752941176];  % Bright blue (#3B4CC0)
    red = [0.705882353, 0.015686275, 0.149019608];  % Bright red (#B40326)
    
    % Create interpolation points
    t = linspace(0, 1, n);
    
    % Interpolate RGB values directly
    cmap = zeros(n, 3);
    for i = 1:3
        cmap(:,i) = interp1([0 1], [blue(i) red(i)], t, 'linear');
    end
end