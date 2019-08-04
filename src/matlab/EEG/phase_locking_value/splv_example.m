% Testing of SPLV on some pre-training lab data (2 trials of 2 individuals)
clear
close
clc

n_chan = 16;
srate = 125;
channels = {'FP1';'FP2';'C3';'C4';'P7';'P8';'O1';'O2';'F7';'F8';'F3';'F4';'T7';'T8';'P3';'P4'};

p1_files = dir('../eeg_data/finger_pointing_2019.06.06/post-training/P16*.csv');
p2_files = dir('../eeg_data/finger_pointing_2019.06.06/post-training/P15*.csv');
chan_loc = '../eeg_data/channel_locs/default_chan_info.ced';

p1_pre1_file = strcat(p1_files(1).folder, '/', p1_files(1).name);
p2_pre1_file = strcat(p2_files(1).folder, '/', p2_files(1).name);

p1_pre1 = eeg_csv_import(p1_pre1_file, chan_loc, 16, srate);
p2_pre1 = eeg_csv_import(p2_pre1_file, chan_loc, 16, srate);

% Before we can run PLV, we need to align the data via their epoch data
p1_epochs = eeg_epochs_from_csv(p1_pre1_file);
p2_epochs = eeg_epochs_from_csv(p2_pre1_file);
data = eeg_time_align(p1_pre1.data, p2_pre1.data, p1_epochs, p2_epochs);

% Attain SPLVs
[results, t] = splv(data, srate, [8 13], 50, 20, 4); 

% Plot normalized SPLVs for electrodes on each subject
figure('Position', [5, 5, 1200, 650])
for i = 1:16
    subplot(4, 4, i);
    sgtitle(sprintf("EEG channel coherence for participants %i and %i in %s", 21,22,"post-training"))
    plot(p1_pre1.times(t)/1000 , squeeze(results(i, i, :)))
    %syncrony_line = refline(0, 1.96);
    %desyncrony_line = refline(0, -1.96);
    %syncrony_line.Color = 'r';
    %desyncrony_line.Color = 'r';
    
    title(sprintf("Coherance in %s", channels{i}))
    xlabel('Time in experiment')
    ylabel('Normalized SPLV')
end

figure('Position', [5, 5, 1200, 650])
for i = 1:16
    subplot(4, 4, i);
    sgtitle(sprintf("EEG channel coherence for participants %i and %i in %s", 21,22,"post-training"))
    plot(p1_pre1.times(t)/1000 , normalize(squeeze(results(i, i, :))))
    syncrony_line = refline(0, 1.96);
    desyncrony_line = refline(0, -1.96);
    syncrony_line.Color = 'r';
    desyncrony_line.Color = 'r';
    
    title(sprintf("Coherance in %s", channels{i}))
    xlabel('Time in experiment')
    ylabel('Normalized SPLV')
end
