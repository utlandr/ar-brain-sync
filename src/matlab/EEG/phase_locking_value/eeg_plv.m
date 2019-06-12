% Here we investigate Praneeth Namburi's implementation of Phase Locking
% Value for an EEG dataset using our own data

% Must run the init.m to setup eeglab first time

clear 
close
clc


% First, import the data
eeg = pop_biosig('../eeg_data/Test-[2019.05.20-11.26.30].gdf');

% specify range of interest
filtSpec.order = 7;
filtSpec.range = [35 45]; %Hz (Gamma range)
lets_go =  pn_eegPLV(eeg, 64, filtSpec);
