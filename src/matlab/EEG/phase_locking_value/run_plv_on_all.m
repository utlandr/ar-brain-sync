% This runs plv on all pairs in a data set. 
%% 0. Init/default settings
n_chan = 16;
srate = 32;
channels = {'FP1';'FP2';'C3';'C4';'P7';'P8';'O1';'O2';'F7';'F8';'F3';'F4';'T7';'T8';'P3';'P4'};
chan_loc = '../eeg_data/channel_locs/default_chan_info.ced';
abs_plv_table = cell2table({"Table", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
abs_plv_table.Properties.VariableNames = ["Pair"; channels];
results_file_name = "vr_first_person_train_abs_results.csv";

%% 1. Define the data pairs to run PLV on
% specify a file containing the file pairs

folder = "../eeg_data/vr_first_person_view/Training/";
pair_file = 'vr_first_person_train.csv';

pairs = readcell(pair_file, 'delimiter' , ',');

%% 2. Loop thorugh each file pairs and run PLV
for i=1:length(pairs)
    try
        p1 = pairs{i,1};
        p2 = pairs{i,2};

        p1_eeg_file = folder+p1;
        p2_eeg_file = folder+p2;

        % Ignore combined files which have a different format
        if contains(p1_eeg_file, "Combined_") || contains(p2_eeg_file, "Combined_")
            continue
        else
            if isfile(p1_eeg_file) && isfile(p2_eeg_file)
                % 2.1 Import actuall eeg data
                p1_eeg_data = eeg_csv_import(p1_eeg_file, chan_loc, 16, srate);
                p2_eeg_data = eeg_csv_import(p2_eeg_file, chan_loc, 16, srate);

                % 2.2 Perform epoch alignment
                p1_epochs = eeg_epochs_from_csv(p1_eeg_file);
                p2_epochs = eeg_epochs_from_csv(p2_eeg_file);
                data = eeg_time_align(p1_eeg_data.data, p2_eeg_data.data, p1_epochs, p2_epochs);
                
                
                % 2.3 Run PLV on data and obtain results
                [results, t] = splv(data, srate, [8 12], 50, 20, 4); 
                
                % 2.4 For each electrode pair, write number of activations
                % PLV above 0.9
                p1_name = strsplit(p1, "_"); 
                group = strsplit(p1_name{end}, ".");
                group = group{1};
                p1_name = p1_name{1};
                p2_name = strsplit(p2, "_");
                p2_name = p2_name{1};
                
                
                pair = strcat(p1_name, "-", p2_name, "-", group);
                ind = size(abs_plv_table,1) +1;
                abs_plv_table.Pair(ind) = pair;
                for j = 1:n_chan
                    abs_res = squeeze(results(j,j,:));
                    abs_count = sum(abs_res > 0.9);
                    
                    abs_plv_table(ind, j+1) = {abs_count};
                end

            else
                warning(sprintf("One or more pair files not found - %s %s", p1_eeg_file, p2_eeg_file));
            end
        end 
    catch exception
        warning(sprintf("Failed PLV run for %s and %s. Error follows:", p1, p2))
        warning(exception.message)
    end

end

%% Output results 
abs_plv_table(1,:) = [];
writetable(abs_plv_table, results_file_name)
