function res = calculate_colour(N, z, c)
    % N is either 1 or 2, corresponding to N1 or N2
    % z is the value of N1 or N2
    % z is either 1, 2 or 3, corresponding to red, green or blue
    
    % scale: from 0 to 200:
    % (-inf, 1], (1, 2], ..., (198, 199], (199, inf]
    %  1          2            199         200
    % categories
    
    % colours: blue if N == 1, red if N == 2:
    % blue: from (r,g,b) = (0.9, 0.9, 1) to (0, 0, 1)
    % red:  from (r,g,b) = (1, 0.9, 0.9) to (1, 0, 0)
    %                       most light      least light
    
    if z <= 1.0
        category = 1;
    elseif z > 199.0
        category = 200;
    else
        category = ceil(z);
    endif
    
    var = 200 - category;
    step = 0.75 / 199.0;
    
    if N == 1
        r = var * step;
        g = var * step;
        b = 1;     
    elseif N == 2
        r = 1;
        g = var * step;
        b = var * step;
    endif
    
    if c == 1
        res = r;
    elseif c == 2
        res = g;
    elseif c == 3
        res = b;
    endif
end