% Testing methods for importing eeg data from csv using eeglab
%
% Author: Reed Bell

% Create options for data import (as original csv is not in table format).
opts = detectImportOptions("pre-training//P1_RW(2019.06.05-11.12.46)__Pre1.csv");
opts.SelectedVariableNames = opts.SelectedVariableNames(3:18);
opts.ExtraColumnsRule = 'ignore';

% View and read matrix.
preview('pre-training/P1_RW(2019.06.05-11.12.46)__Pre1.csv', opts)
eeg_data_test = readmatrix('pre-training/P1_RW(2019.06.05-11.12.46)__Pre1.csv', opts);