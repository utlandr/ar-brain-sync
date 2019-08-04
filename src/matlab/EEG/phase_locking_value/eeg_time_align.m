function [eeg_data, epoch_range] = eeg_time_align(eeg1, eeg2, eeg1_epochs, eeg2_epochs)
    if isequal(eeg1_epochs, eeg2_epochs)
        eeg_data(:, :, 1) = eeg1(:, :);
        eeg_data(:, :, 2) = eeg2(:, :);
        epoch_range = minmax(eeg1_epochs);
        
    else
        eeg1_epoch_set = unique(eeg1_epochs);
        eeg2_epoch_set = unique(eeg2_epochs);

        fin_ind = find(ismember(eeg1_epoch_set,eeg2_epoch_set), 2, 'last');
        start_ind = find(ismember(eeg1_epoch_set,eeg2_epoch_set), 2);

        if isempty(fin_ind) || isempty(start_ind)
            error("eeg datasets have no overlapping epochs")

        else
            start_epoch = eeg1_epoch_set(start_ind(2));
            fin_epoch = eeg1_epoch_set(fin_ind(1));
        
            eeg1_data_ind = find(eeg1_epochs == start_epoch, 1):find(eeg1_epochs == fin_epoch, 1, 'last');
            eeg2_data_ind = find(eeg2_epochs == start_epoch, 1):find(eeg2_epochs == fin_epoch, 1, 'last');
            if length(eeg1_data_ind) ~= length(eeg2_data_ind)
                error("Dataset sizes after time adjust do not match")
    
            else
                eeg_data(:, :, 1) = eeg1(:, eeg1_data_ind);
                eeg_data(:, :, 2) = eeg2(:, eeg2_data_ind);
                epoch_range = [start_epoch fin_epoch];
            
        end
    end
    
    
    
end



