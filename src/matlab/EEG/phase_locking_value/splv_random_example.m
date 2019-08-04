% Testing of SPLV on randomly generated data
clear 
close
clc

for i = 1:10
    
    % Random data
    eeg_rand_1 = rand(16,1250,10);

    % options
    srate = 125;
    freq_band = [8 12];
    order = 50;
    delta_half = 10;
    
    % Calculate plv data
    splv_output = splv(eeg_rand_1, srate, freq_band, order, delta_half, 19);

    eeg_all(i, :) = squeeze(splv_output(1, 1, :));
end

% Plot comparison of node 1
plot(mean(eeg_all));