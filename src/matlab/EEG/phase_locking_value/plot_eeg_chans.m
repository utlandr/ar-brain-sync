function plot_eeg_chans(eeg1, eeg2, n_channels)

    % Plot data at an electrode for over time
    figure('Position', [5, 5, 1200, 650])
    for i = 1:n_channels
        subplot(4, 4, i);
        series1 = eeg1.data(i, :);
        series2 = eeg2.data(i, :);
        plot(eeg1.times, series1);
        hold on;
        plot(eeg2.times ,series2);
        hold off;

        xlabel('Time')
        ylabel('EGG reading (uV)')

        title(i)
    end

end