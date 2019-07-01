function [plv, t_points] = splv(eeg_data, srate, freq_band, order, delta_half, overlap)
    % DOCSTRING
    % This is an implementation of the Phase Locking Value (PLV) algorithm that
    % only considers a single sample/trial to compute PLVs. It varies from the 
    % traditional implementation which averages over the samples, by averaging
    % over a time window. This implmentation is similar to pn_eegPLV
    %
    % INPUT
    % eeg_data      -   (nxmx2 double)A 3D matrix with dimensions nchannels 
    %                   x ntimepoints x ntrials = 2 of EEG readings.
    % srate         -   (int) The sampling rate of the EGG data in Hz
    % freq_band     -   (2x1 double) The range of frequencies to calculate PLV
    %                   for (e.g. use [8 12] for alpha waves.
    % order         -   (int) The order of the FIR filter (should be designed to
    %                   approximately 4-5 cycles of the desired frequency. e.g
    %                   alpha wave has frequency of 10Hz (100ms) with a srate
    %                   of 125Hz (8ms) would use an order = 50 so that 8*50 =
    %                   400ms/100ms = 4 cycles
    % overlap       -   (int) specifies the number of data points of overlap
    %                   at each end of the of the time window.
    %
    % OUTPUT
    % plv           -   (nxmx2 double) A 3D matrix with dimensions nchannels x
    %                   ntimepoints x ntrials = 2 containing calculated PLVs
    %
    % AUTHOR
    % Reed Bell     -   rbel068@aucklanduni.ac.nz
    
    % First, recieve information about channel number and timepoints
    [n_channels, n_timepoints, ~] = size(eeg_data);
    
    % Specify FIL filter for specific frequency range
    fil = fir1(order, 2/srate*freq_band);
                    
    % Apply the FIL filter
    fil_eeg = filter(fil, 1, eeg_data, [], 2);
    fil_hilbert_eeg = fil_eeg;
    
    % Apply Hilbert transformation to find instantaneous phase at each
    % timepoint
    for k = 1:2 
        for channel = 1:n_channels
            tmp_hilbert = hilbert(fil_eeg(channel, :, k));
            fil_hilbert_eeg(channel, :, k) = angle(tmp_hilbert);
            
        end
    end
    
    % Calculate time averaged PLV between all possible channel pairs
    t_points = 1:2*delta_half-overlap:n_timepoints;
    plv = zeros(n_channels, n_channels, length(t_points));
    
    for chan_1 = 1:n_channels 
        for chan_2 = 1:n_channels
            plv_ind = 1;
            for t = t_points               
                % Setting the time window
                if t - delta_half + 1 < 1
                    window = 1:t+delta_half;
                    
                elseif t + delta_half > n_timepoints
                    window = t-delta_half+1:n_timepoints;
                    
                else
                    window = t-delta_half+1:t+delta_half;
                    
                end
                
                single_plvs = exp(1i*(fil_hilbert_eeg(chan_1, window, 1) - fil_hilbert_eeg(chan_2, window, 2)));
                plv(chan_1, chan_2, plv_ind) = abs((1 / (2 * delta_half))*sum(single_plvs));
                plv_ind = plv_ind+1;
            
            end
        end
    end
    
    
    return