% Testing of SPLV on some pre-training lab data (2 trials of 2 individuals)
clear
close
clc

p1_files = dir('../eeg_data/finger_pointing_2019.06.06/post-training/P21*.csv');
p2_files = dir('../eeg_data/finger_pointing_2019.06.06/post-training/P22*.csv');
chan_loc = '../eeg_data/channel_locs/default_chan_info.ced';
n_chan = 16;
srate = 125;

p1_pre1_file = strcat(p1_files(1).folder, '/', p1_files(1).name);
p2_pre1_file = strcat(p2_files(1).folder, '/', p2_files(1).name);

p1_pre1 = eeg_csv_import(p1_pre1_file, chan_loc, 16, srate);
p2_pre1 = eeg_csv_import(p2_pre1_file, chan_loc, 16, srate);

% Before we can run PLV, we need to align the data via their epoch data
p1_pre1_epochs = eeg_epochs_from_csv(p1_pre1_file);
p2_pre1_epochs = eeg_epochs_from_csv(p2_pre1_file);



data = zeros(n_chan, fin_p1 - start_p1 + 1, 2);
data(:, :, 1) = p1_pre1.data(:, start_p1:fin_p1);
data(:, :, 2) = p2_pre1.data(:, start_p2:fin_p2);

[results, t] = splv(data, srate, [20 30], 50, 10, 2);

    % Plot data at an electrode for over time
    figure('Position', [5, 5, 1200, 650])
    for i = 1:16
        subplot(4, 4, i);
        plot(t , squeeze(results(i, i, :)))

        xlabel('Time')
        ylabel('SPLV')

        title(i)
    end

