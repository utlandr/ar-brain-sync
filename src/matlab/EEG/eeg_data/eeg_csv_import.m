
function eeg_dataset = eeg_csv_import(csv_file, chan_file, nbchan, srate)
    % DOCSTRING
    % Import a csv file containign eeg information (from inlab experiment) and
    % convert into a eeglab compliat dataset for data analysis
    %
    % INPUT
    % csv_file      -   (char array) the relative location of csv file
    % chan_file     -   (char array) the relative location of  a supported 
    %                   channel file
    % nbchan        -   (int) the number of channels in the dataset
    % srate         -   (int) the recording frequency of the headset (e.g. 125)
    %
    % OUTPUT   
    % eeg_dataset   -   (struct) an eeglab compliant dataset
    % 
    % AUTHORS 
    % Reed Bell     -   rbel068@aucklanduni.ac.nz    
    % Gus Stone

    % Create options for data import (as original csv is not in table format). 
    opts = detectImportOptions(csv_file);
    opts.ExtraColumnsRule = 'ignore';
    eeg_data_opts = opts;
    eeg_time_epoch_opts = opts;

    % Setup options for eeg data
    eeg_data_opts.SelectedVariableNames = opts.VariableNames(3:18);

    % Options for time/epoch information
    eeg_time_epoch_opts.SelectedVariableNames = opts.VariableNames(1:2);

    % View and read matrix (trnaspose for eeglab import).
    eeg_data = transpose(readmatrix(csv_file, eeg_data_opts));
    eeg_time_epoch = transpose(readmatrix(csv_file, eeg_time_epoch_opts));
    
    start_time = eeg_time_epoch(1,1);
    
    % Place this data into an eeglab struct containing all relevant info
    eeg_dataset = pop_importdata('setname', '', ...
                                 'data', eeg_data, ...
                                 'dataformat', 'array', ...
                                 'nbchan', nbchan, ...
                                 'chanlocs', chan_file, ...
                                 'srate', srate, ...
                                 'xmin', start_time);
    
end