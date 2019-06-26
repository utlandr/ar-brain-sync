% Here we investigate Praneeth Namburi's implementation of Phase Locking
% Value for an EEG dataset using our own data

% Must run the init.m to setup eeglab first time

clear 
close
clc


% First, import the data
eeg1 = pop_biosig('../eeg_data/p3.gdf');
egg2 = pop_biosig('../eeg_data/reference/p4.gdf');

% specify range of interest
filtSpec.order = 7;
filtSpec.range = [10 30]; %Hz
plv_ouput =  pn_eegPLV(eeg.data, 64, filtSpec);
