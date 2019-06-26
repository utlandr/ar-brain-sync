close all

% Running the example PLV given in the pn_eegPLV example
eegData = rand(16, 125*10, 10); 
srate = 125; %Hz
filtSpec.order = 1;
filtSpec.range = [35 45]; %Hz
dataSelectArr = rand(10, 1) >= 0.5; % attend trials
dataSelectArr(:, 2) = ~dataSelectArr(:, 1); % ignore trials
plv = pn_eegPLV(eegData, srate, filtSpec, dataSelectArr);
figure; plot((0:size(eegData, 2)-1)/srate, squeeze(plv(:, 1, 2, :)));
xlabel('Time (s)'); ylabel('Plase Locking Value');

dataSelectArr