% Here, we attempt to apply PCA to EEG data to identify activations in the
% data. For now, we will rely on a single dataset, but we can look at more
% in the future.


close all 
clear all
clc

% First, run `eeglab` to import eeglab functions

% Import EEG data
EEG = pop_biosig('p3.gdf');
[n_var, n_epoch] = size(EEG.data);

c={ 'FP1','FP2','C3','C4', 'P7','P8','O1','O2', 'F7','F8','F3','F4', 'T7','T8','P3','P4' };

window_size = 5;

% create a time window in which we apply PCA
for i = 400:600%1:5:n_epoch - window_size
    eeg_tmp = EEG.data(:, i:i+5);
    pca_tmp = pca(eeg_tmp');
    
    figure(1)
    scatter(pca_tmp(:, 1), pca_tmp(:, 2));
    
    figure(2)
    scatter(eeg_tmp(:, 1), egg_tmp(:, 2))
    
    pause(2)
    
    
    
    
end

