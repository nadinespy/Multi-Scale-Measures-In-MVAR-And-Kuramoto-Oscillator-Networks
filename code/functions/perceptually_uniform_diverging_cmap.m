function cmap = perceptually_uniform_diverging_cmap(n)
    % Creates a perceptually uniform diverging colormap with enhanced contrast
    % n: number of color steps (should be odd for symmetric colormap)
    
    if nargin < 1
        n = 255;
    end
    % Ensure n is odd for symmetric colormap
    if mod(n, 2) == 0
        n = n + 1;
    end
    
    % Define endpoints in LAB space with enhanced contrast
    lab_neg = [30 85 -85];  % More saturated blue
    lab_mid = [95 0 0];     % Brighter white
    lab_pos = [30 85 45];   % More saturated red
    
    % Create interpolation points
    half_n = floor(n/2);
    t_neg = linspace(0, 1, half_n + 1);
    t_pos = linspace(0, 1, half_n + 1);
    
    % Interpolate in LAB space
    lab_lower = zeros(half_n + 1, 3);
    lab_upper = zeros(half_n + 1, 3);
    for i = 1:3
        lab_lower(:,i) = interp1([0 1], [lab_neg(i) lab_mid(i)], t_neg, 'pchip');
        lab_upper(:,i) = interp1([0 1], [lab_mid(i) lab_pos(i)], t_pos, 'pchip');
    end
    
    % Combine lower and upper parts
    lab_combined = [lab_lower(1:end-1,:); lab_upper];
    
    % Convert from LAB to RGB (same conversion code as sequential function)
    cmap = zeros(size(lab_combined));
    for i = 1:size(lab_combined, 1)
        % First convert LAB to XYZ
        L = lab_combined(i,1);
        a = lab_combined(i,2);
        b = lab_combined(i,3);
        
        % LAB to XYZ conversion
        fy = (L + 16) / 116;
        fx = a / 500 + fy;
        fz = fy - b / 200;
        
        % Reference white point (D65)
        Xn = 0.95047;
        Yn = 1.00000;
        Zn = 1.08883;
        
        % Convert fx, fy, fz to X, Y, Z
        if fy > 6/29
            Y = Yn * fy^3;
        else
            Y = (fy - 16/116) * 3 * (6/29)^2 * Yn;
        end
        if fx > 6/29
            X = Xn * fx^3;
        else
            X = (fx - 16/116) * 3 * (6/29)^2 * Xn;
        end
        if fz > 6/29
            Z = Zn * fz^3;
        else
            Z = (fz - 16/116) * 3 * (6/29)^2 * Zn;
        end
        
        % XYZ to RGB conversion matrix (sRGB)
        M = [3.2406 -1.5372 -0.4986;
            -0.9689 1.8758 0.0415;
            0.0557 -0.2040 1.0570];
        rgb = M * [X; Y; Z];
        
        % Apply gamma correction and clip values
        mask = rgb <= 0.0031308;
        rgb(mask) = 12.92 * rgb(mask);
        rgb(~mask) = 1.055 * rgb(~mask).^(1/2.4) - 0.055;
        
        % Clip values to [0,1] range
        rgb = min(max(rgb, 0), 1);
        cmap(i,:) = rgb';
    end
end
% Example usage:
% n = 255;  % number of colors
% cmap = perceptualDiverging(n);
% colormap(cmap);
% imagesc(your_data);
% colorbar;

% Optional: Add this to ensure proper scaling of the colorbar
% clim = max(abs(your_data(:)));
% set(gca, 'CLim', [-clim clim]);
