tpoints = 1:50;
d_half = 10;
n_overlap = 2;
n_timepoints = length(tpoints);

for t = 1:2*d_half-n_overlap:n_timepoints
    if t - d_half + 1 < 1
        window = 1:t+d_half;
                    
    elseif t + d_half > n_timepoints
        window = t-d_half+1:n_timepoints;

    else
        window = t-d_half+1:t+d_half;

    end
    
    window
end