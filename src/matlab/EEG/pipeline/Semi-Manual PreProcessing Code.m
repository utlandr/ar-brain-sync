%% NOTES

% At this stage the files are already cut up (we need to automate that
% prior process. Data will be in EDF format as well.

% Also need plugin for hdf5 format on EEGLAB (as files may originally be in
% hdf

% SET/UNSET the dataset options for generation the dataset. Set will give
% you just one. The other will give two (.set and .fdt) 

% eegh is a godsend. will show all the things you do in gui

%%UPDATE 29/04/19
% Gamma di ganti freq 40
%% Isilah titik titik dibawah ini
run E:\eeglab14_1_2b\eeglab; %EEG home
nchannels = 14 ; % num of channels in the EEG headset
ChanLocs = 'E:\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp'; %STANDARD chan locs
ChanLocs2 = 'E:\eeglab14_1_2b\emotivupdated.ced'; % path for custom channel locs
ChanNum = [1:16]; % starting from 1
% ChanNum = [1:19]; %if Deymed / 19 Channels Stuffs
CutPath = 'C:\Users\DAI 02 - Neurolab\Desktop\Cutting Folder\CutFiles\'; % placing all the cut files in this folder
close;

%STEP 0: Change the option to use double precision
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1,...
        'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0,...
        'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1,...
        'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
        'option_checkversion', 1, 'option_chat', 0);

%% PREPROCESSING
for SubjID = 1: length(SetFiles);
    loadName = SetFiles(SubjID).name;
    dataName = loadName(1:end-4);
    saveName = [dataName '_PP'] %adding 'PP' (PreProcessed) to the end 
    
    %Step 1: Import data. REVISED 
  % EEG = pop_biosig(loadName, CutPath); %IF .EDF % loading hd5
    EEG = pop_loadset(loadName, CutPath); %IF .SET
    EEG.setname = dataName;

    %STEP 2: Changing sample rate
    EEG = pop_resample(EEG, 500); % changed to 500

    %STEP 3: High-pass filter and Low pass filter
    EEG = pop_eegfiltnew(EEG,1,45,826); % What does 826 do

    %STEP 4: Select channels
    EEG = pop_select( EEG,'channel',ChanNum);
    
    %STEP 5: Importing channel location
    EEG = pop_chanedit(EEG, 'lookup',ChanLocs,'load',{ChanLocs2 'filetype' 'autodetect'}); % What does second loc file do
    
        %Keeping original EEG file
        originalEEG = EEG; % Just back up original EEG data
        
    %STEP 6: Removing bad channels by using raw_clean plug in
    EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5); % eeg is dirty we clean it here using raw_clean plugin

    %STEP 7: Interpolate channels.
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical'); %we interpolate bad data from above to replace it (and maintain electrode data).

    %STEP 8: Apply average reference after adding initial reference
    EEG.nbchan = EEG.nbchan+1; % Important in all electrodes. Must find the average to use as reference
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []); % pop_reref does the averaging, where does it go?
    EEG = pop_select( EEG,'nochannel',{'initialReference'}); % back to the start

    %STEP 9: Epoching data 1 to 3 sec 
    EEG = eeg_regepochs(EEG, 'limits', [1 2] , 'extractepochs', 'on'); % creating epochs for FFT (NOT PLV)... What do the options do?
    
    %STEP 10: Automatic epoch rejection
    EEG = pop_autorej(EEG, 'threshold', 1000,'startprob',5,'maxrej', 5, 'nogui','on'); % Automaitcally finds the noise and removes the epochs that suck
    
    %STEP 11: Rejection epoch by probability (6SD single channel, 2SD for all channels)
    % manual rejection if data in each chan is 6StanDevs away.
    EEG = eeg_checkset( EEG ); % what happens after removal of epochs, what is joint prob
    EEG = pop_jointprob(EEG,1,[1:16] ,6,2,1,0,0,[],0); % Why manual removal
    
    %STEP 12: Running ICA
    % seperates the contributions at each electrode.
    EEG = eeg_checkset( EEG );
    EEG = pop_runica(EEG, 'extended',1,'interupt','on'); % what are extended and interupt do?
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    %STEP 13 : Checking whether EEG data contains ICA decomposition
    EEG = eeg_checkset(EEG, 'ica');
    
    %STEP 14: Removing line noise with cleanLine;
    % electricity interference etc (requires plugin)
    EEG = pop_cleanline(EEG,'Bandwidth',2,'ChanCompIndices',[1:EEG.nbchan],...
        'ComputeSpectralPower',0,'LineFrequencies',[60 120],...
        'NormalizeSpectrum',0,'LineAlpha',0.01,'PaddingFactor',2,...
        'PlotFigures',0,'ScanForLines',1,'SignalType','Channels',...
        'SmoothingFactor',100,'VerboseOutput',1); % HAPUS SlidingWinLength dan SlidingWinstep
    
    %STEP 15: Rejecting ICA by extreme value
    %detects artifacting (and rejects the ICA)
    EEG = pop_eegthresh(EEG,0,[1:5] ,-20,20,1,2.996,0,1); % double check start stop times

    %STEP 16: Finding Power for each frequency
    %%%%%%%%%%%%%%%%%%% Finding power for each frequency band %%%%%%%%%%%%%%%%%%%
    
    EEG = pop_saveset(EEG, 'filename',saveName,'filepath', ProcPath);  
end

%% Calculate the PLV
% Does PLV use components, or the channel data for PLV.

%% Statistic
for ProcID = 1:length(ProcFiles)
    loadProc = ProcFiles(ProcID).name;
    procData = loadProc(1:end-4);
    
    EEG = pop_loadset(loadProc, ProcPath); %IF .SET
    EEG.setname = procData;
    
    for n = 1: nchannels

        [spectra,freqs] = spectopo(EEG.data(n,:,:), 0, EEG.srate, 'plot', 'off'); % Sesuaikan channel mana yang mau diambil

        % delta=1-4, theta=4-8, alpha=8-13, beta=13-30, gamma=30-80
        deltaIdx{n} = find(freqs>1 & freqs<=4);
        thetaIdx{n} = find(freqs>4 & freqs<=8);
        alphaIdx{n} = find(freqs>8 & freqs<=13);
        betaIdx{n}  = find(freqs>13 & freqs<=30);
        gammaIdx{n} = find(freqs>30 & freqs<=80);

        % compute absolute power
        deltaPower{n} = mean(10.^(spectra(deltaIdx{n})/10));
        thetaPower{n} = mean(10.^(spectra(thetaIdx{n})/10));
        alphaPower{n} = mean(10.^(spectra(alphaIdx{n})/10));
        betaPower{n}  = mean(10.^(spectra(betaIdx{n})/10));
        gammaPower{n} = mean(10.^(spectra(gammaIdx{n})/10));
    end
    
        deltaCell= cell2mat(deltaPower);
        thetaCell= cell2mat(thetaPower);
        alphaCell= cell2mat(alphaPower);
        betaCell = cell2mat(betaPower);
        gammaCell= cell2mat(gammaPower);
        
        % Dikumpul per wave
        deltaAll(ProcID,:) = deltaCell
        thetaAll(ProcID,:) = thetaCell
        alphaAll(ProcID,:) = alphaCell
        betaAll (ProcID,:) = betaCell
        gammaAll(ProcID,:) = gammaCell
        
end
%Indexing p-value per Domain
for StatID = 1:JumlahDomain+1:length(ProcFiles);
%     loadStat = ProcFiles(StatID).name;
%     statData = loadStat(1:end-4);
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(deltaAll(StatID,:),deltaAll(StatID+Stat,:));
     PValDelta(StatID+Stat) = p ;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(thetaAll(StatID,:),thetaAll(StatID+Stat,:));
     PValTheta(StatID+Stat) = p;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(alphaAll(StatID,:),alphaAll(StatID+Stat,:));
     PValAlpha(StatID+Stat) = p;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(betaAll(StatID,:),betaAll(StatID+Stat,:));
     PValBeta(StatID+Stat) = p;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(gammaAll(StatID,:),gammaAll(StatID+Stat,:));
     PValGamma(StatID+Stat) = p;
     end
end
     %removing zeroes for Indexing
PValDelta(PValDelta==0) =[];
PValTheta(PValTheta==0) =[];
PValAlpha(PValAlpha==0) =[];
PValBeta(PValBeta==0) =[];
PValGamma(PValGamma==0) =[];

%Dikumpul jadi satu workspace also for 1 row Indexing
PValSubj = [PValDelta; PValTheta; PValAlpha; PValBeta; PValGamma]
[PValMin,PValIdx] = min(PValSubj)

%% Plotting
for z = 1:length(PValIdx);
    loadStat = ProcFiles(z).name;
    statData = loadStat(1:end-4);
 
%indexing smallest p-value sebagai domain (1=delta - 5=gamma)
 if PValIdx(z) == 1
           Coord(z) = {deltaFrame}
           ImgName{z} = ['_delta.png']
        elseif PValIdx(z) == 2
           Coord(z) =  {thetaFrame}
           ImgName{z} = ['_theta.png']
        elseif PValIdx(z) == 3
           Coord(z) =  {alphaFrame}
           ImgName{z} = ['_alpha.png']
        elseif PValIdx(z) == 4
          Coord(z) =  { betaFrame}
          ImgName{z} = ['_beta.png']
  elseif PValIdx(z) == 5
          Coord(z) = {gammaFrame}
          ImgName{z} = ['_gamma.png']
          
 end
end
%move BaseFile to new Folder
for BaseSet = 1:JumlahDomain+1:length(ProcFiles)
movefile(ProcFiles(BaseSet).name, ProcBasePath)
movefile (ProcFDT(BaseSet).name, ProcBasePath)
end
ProcFiles = dir([ProcPath '*.set']); %refresh filelist in the folder after moved
  for PlotID = 1:length(ProcFiles)
      PlotFile = ProcFiles(PlotID).name
      PlotName = [PlotFile(1:end-4),ImgName{PlotID}]
            EEG = pop_loadset(PlotFile, ProcPath);
      EEG.setname = PlotName;
        %open figure&
        figure('units','normalized','outerposition',[0 0 1 1]);  %%Full Screen
        title(PlotName); %Title
        % Plotting Head
        pop_spectopo(EEG, 1, [], 'EEG', ...
        'percent', 100, 'freq', [2 6 12 25 40], 'freqrange',[0 45], 'electrodes','off', 'overlap', 0);
        set(gcf,'color','w');
        %Theta
        F = getframe(gcf, Coord{PlotID});
        figure;
        imshow(F.cdata);
        Image = frame2im(F);
               imwrite(Image,[ImgPath PlotName],'png')
        close
        close
   end