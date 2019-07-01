function epochs = eeg_epochs_from_csv(csv_file)
    opts = detectImportOptions(csv_file);
    opts.ExtraColumnsRule = 'ignore';
    opts.SelectedVariableNames = 'Epoch';
    epochs = readmatrix(csv_file, opts);
    

end