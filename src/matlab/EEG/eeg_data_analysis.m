% eeg_data_analysis.m
%
% Description: A quick look into what the EEG data from a finger tracing
% excercise looks like.
%
% Author: Reed Bell

close all
clc

% Import data
EEG1=pop_biosig('p3.gdf');
EEG2=pop_biosig('p4.gdf');

% Channel names
chan={ 'FP1', 'FP2', 'C3', 'C4', ...
       'P7', 'P8', 'O1', 'O2', ...
       'F7', 'F8', 'F3', 'F4', ...
       'T7', 'T8', 'P3', 'P4' };

% Generate distribution data for a single electrode on both datasets
width = 0.5;
figure('Position', [5, 5, 1200, 650])

for i = 1:length(chan)
    subplot(4, 4, i);
    series1 = EEG1.data(i, :);
    series2 = EEG2.data(i, :);
    hist1 = histogram( series1,  'EdgeColor', 'blue', 'FaceColor',  'blue', 'BinWidth', width);
    hold on;
    hist2 = histogram( series2, 'EdgeColor', 'green', 'FaceColor',  'green', 'FaceAlpha', 0.2 , 'BinWidth', width);
    hold off;
    
    xlabel('EEG Average Reading')
    ylabel('Frequency')
    
    title(chan{i})
end

% NOTE: Results display that bot datasets have completely different ranges:
%       - EEG1 (p3.gdf) displays majority readings between 14.8 and 15.1.
%         Some readings (notably P4) display low/near zero readings which
%         is odd.
%       
%       - EEG2 (p4.gdf) displays majority readings between 5-10. Notably,
%         F8 is an exception here where readings are consistantly above 10.


figure(2)
histogram(EEG1.data(12, :), 'BinWidth', 0.001)
xlabel('Average EEG Reading')
ylabel('Frequency')
title('EEG1 data distribution at electrode F4')
figure(3)
histogram(EEG2.data(12, :), 'BinWidth', 0.1)
xlabel('Average EEG Reading')
ylabel('Frequency')
title('EEG2 data distribution at electrode F4')

% Plot data at an electrode for over time
figure('Position', [5, 5, 1200, 650])
for i = 1:length(chan)
    subplot(4, 4, i);
    series1 = EEG1.data(i, :);
    series2 = EEG2.data(i, :);
    line1 = plot(1:EEG1.pnts, series1);
    hold on;
    line2 = plot(1:EEG2.pnts ,series2);
    hold off;
    
    xlabel('Epoch')
    ylabel('EGG average reading')
    
    title(chan{i})
end




