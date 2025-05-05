function newmap = bluered(m)
    % Define the three main colors with more intense blue and red
    blue = [0 0 1];      % Pure saturated blue
    white = [1 1 1];     % White middle
    red = [1 0 0];       % Pure saturated red
    
    % Create simple interpolation
    new = [blue; white; red];
    oldsteps = linspace(0, 1, 3);  % Linear spacing for smoother transition
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);

    %new = [blue; blue; white; red; red];  % Duplicate end points for sharper transition
    %oldsteps = [0 0.5 1 1.5 2];         % Adjust spacing to make transitions sharper
    %newsteps = linspace(0, 1, m);
    %newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
end