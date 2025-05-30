

datapath = ('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/');
addpath('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/seeg'); % dont genpath
addpath('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/resources/fieldtrip-20240110'); % dont genpath
addpath(('/Users/shanemckeon/Documents/GitHub/eeglab')) % dont genpath
addpath '/Users/shanemckeon/Downloads/uANTS_synch_MATLABfiles'
addpath(genpath('/Users/shanemckeon/Documents/GitHub/generalCode'))
addpath(genpath('/Users/shanemckeon/Documents/GitHub/sEEG/usefulCode'))


setfiles0 = dir([datapath, 'P*/rest/rest_eeglabFormat_referenced.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

% for i = 1:length(setfiles)
i=14;
%% clean seeg electodes
pt = (setfiles0(i).folder(106:112));
inputfile = (fullfile([setfiles0(i).folder,'/', setfiles0(i).name]));

[badChans] = findNoisyChans(pt, 'rest'); 

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG = pop_loadset(inputfile); % load in eeg file

% Specify bad channels to remove
bad_channels = badChans{1,1};

% Find indices of bad channels
bad_channel_indices = find(ismember({EEG.chanlocs.labels}, bad_channels));

% Remove bad channels
EEG = pop_select(EEG, 'nochannel', bad_channel_indices);

%% preproc scalp if they have them

EEG = preprocessScalpEEG(EEG);


%end