
 addpath /Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/usefulCode
 addpath /Users/shanemckeon/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/Colorspace Transformations
 addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/generalCode'))

% addpath(genpath('/Users/shanemckeon/Documents/GitHub/generalCode'))
% addpath(genpath('/Users/shanemckeon/Documents/GitHub/sEEG/usefulCode'))

addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/resources/eeglab2024.0'));
eeglab

addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/resources/fieldtrip-20240110'))
ft_defaults

datapath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/';
outpath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/';

% resultPath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/individual_subject_files/';
% photoPath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/channelFigures/';

resultPath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/acw/individual_subject_files/';
photoPath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/channelFigures/';


setfiles0 = dir([datapath, 'P*/rest/rest_eeglabFormat_referenced.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

%% plot for every person individually

for i = 1:length(setfiles0)

    pt = (setfiles0(i).folder(106:112));

    if ~isfile([photoPath, pt, '.png'])

        inputfile = (fullfile([setfiles0(i).folder,'/', setfiles0(i).name]));

        % open eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

        EEG = pop_loadset(inputfile); % load in eeg file

        channels = num2cell(zeros(length(EEG.chanlocs), 1));
        channel_labels = zeros(length(EEG.chanlocs), 1);

        for c = 1:length(EEG.chanlocs)
            channels{c,:} = {EEG.chanlocs(c).labels};
        end

        channel_labels = cellfun(@(x) x{1}, channels, 'UniformOutput', false);

        % Find rows where labels start with '_'
        rows_to_remove = startsWith(channel_labels, '_');

        % Remove those rows
        channel_labels(rows_to_remove) = [];

        plotInflSurf('pt',pt, channel_labels);

        saveas(gcf, [photoPath, pt, '.png']);
        close all;

    end
end

%% All channels, one plot

subject_ids = cell(length(setfiles0), 1);
all_channels = cell(length(setfiles0), 1);

for i = 1:length(setfiles0)

    pt = setfiles0(i).folder(106:112);  % Extract subject ID

    inputfile = fullfile(setfiles0(i).folder, setfiles0(i).name);

    % Open EEGLAB
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    % Load EEG file
    EEG = pop_loadset(inputfile);

    % Extract channel labels
    channel_labels = cell(length(EEG.chanlocs), 1);

    for c = 1:length(EEG.chanlocs)
        channel_labels{c} = EEG.chanlocs(c).labels;
    end

    % Remove channel labels that start with '_'
    %channel_labels(startsWith(channel_labels, '_')) = [];

    % Store subject ID and channel labels
    subject_ids{i} = pt;
    all_channels{i} = channel_labels;

end

%% color per person
num_patients = length(setfiles0);
patient_colors = linspace(-1, 1, num_patients); % Assigns a unique scalar value per patient

% Initialize the cell array for channel colors
channel_colors = cell(size(all_channels));

% Assign a unique scalar color to each patient's channels
for i = 1:num_patients
    num_channels = length(all_channels{i});
    channel_colors{i} = repmat(patient_colors(i), num_channels, 1); % Assign the same color to all channels of a patient
end


% Call the function with properly formatted input
plotInflSurf('mni', subject_ids, all_channels, channel_colors);

%% color by lobe 

lobeTable = readtable('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/allSubChannelsAssignedToLobe.csv');

uniqueSubs = unique(lobeTable.Subject);

% Initialize cell arrays
subjectArray = cell(1, numel(uniqueSubs));
chanLabelArray = cell(1, numel(uniqueSubs));

% Loop through each subject and extract corresponding Channel_Label values
for i = 1:numel(uniqueSubs)
    thisSub = uniqueSubs{i};
    idx = strcmp(lobeTable.Subject, thisSub);
    
    subjectArray{i} = thisSub;  % Store subject ID
    chanLabelArray{i} = lobeTable.Channel_Label(idx);  % Store channel labels for that subject
end

% Initialize cell array to hold lobe codes per subject
lobeCodeArray = cell(1, numel(uniqueSubs));

% Loop through each subject
for i = 1:numel(uniqueSubs)
    thisSub = uniqueSubs{i};
    idx = strcmp(lobeTable.Subject, thisSub);
    
    % Get lobe labels for this subject
    lobes = lobeTable.Lobe(idx);
    
    % Convert lobe labels to numeric codes
    codes = zeros(size(lobes,1),1);
    for j = 1:numel(lobes)
        switch lower(strtrim(lobes{j}))
            case 'frontal'
                codes(j,:) = 1;
            case 'parietal'
                codes(j,:) = 2;
            case 'occipital'
                codes(j,:) = 3;
            otherwise
                codes(j,:) = NaN; % Optional: handle unexpected labels
        end
    end
    
    lobeCodeArray{i} = codes;
end

% Call plotInflSurf once with all channels and their colors
plotInflSurf('mni', subjectArray, chanLabelArray, lobeCodeArray,[],[],[],[0.0039, 0.1686, 0.0510; 0.1529, 0.4627, 0.2392; 0.5647, 0.7216, 0.6510]);

saveas(gcf, [photoPath, 'lobeColors_left.pdf']);




%% color by acw rank

% acw = readtable('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/allSubjects_ACW_avgEpochs.csv');
acw = readtable('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/acw/allSubjects_ACW_avgEpochs.csv');

acwValues = table2array(acw(:,4));

ranks = tiedrank(acwValues, 'descend'); % rank the acw values

rankValue = [ranks acwValues];

% divide the ranks for each person
lastRow = 1;

for s = 1:length(setfiles0)
    numChannels = size(channel_colors{s},1)-1;
    subACWrank{s} = ranks((lastRow:lastRow+numChannels), :);

    lastRow = lastRow + numChannels;

end

plotInflSurf('mni', subject_ids, all_channels, subACWrank);


%% color by acw rank adults only

acw = readtable('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/allSubjects_ACW_avgEpochs.csv');
acw_filtered = acw(acw.age >= 18, :);

filtered_subject_ids = acw_filtered.Subject;


subject_ids = cell(length(setfiles0), 1);
all_channels = cell(length(setfiles0), 1);

for i = 1:length(setfiles0)

    pt = setfiles0(i).folder(106:112);  % Extract subject ID
    if ismember(pt, filtered_subject_ids)

        inputfile = fullfile(setfiles0(i).folder, setfiles0(i).name);

        % Open EEGLAB
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

        % Load EEG file
        EEG = pop_loadset(inputfile);

        % Extract channel labels
        channel_labels = cell(length(EEG.chanlocs), 1);

        for c = 1:length(EEG.chanlocs)
            channel_labels{c} = EEG.chanlocs(c).labels;
        end

        % Remove channel labels that start with '_'
        channel_labels(startsWith(channel_labels, '_')) = [];

        % Store subject ID and channel labels
        subject_ids{i} = pt;
        all_channels{i} = channel_labels;
    end

end

non_empty_indices = find(~cellfun('isempty', subject_ids));

subject_ids_adults = subject_ids(~cellfun('isempty', subject_ids));
all_channels_adults = all_channels(~cellfun('isempty', all_channels));

subACWrank_adults = subACWrank(non_empty_indices);


plotInflSurf('mni', subject_ids_adults, all_channels_adults, subACWrank_adults);

saveas(gcf, [photoPath, 'allChannels_ACWColor_AdultsOnly_left.pdf']);



%% color by acw rank adolescence only

acw = readtable('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/allSubjects_ACW_avgEpochs.csv');
acw_filtered = acw(acw.age < 18, :);

filtered_subject_ids = acw_filtered.Subject;


subject_ids = cell(length(setfiles0), 1);
all_channels = cell(length(setfiles0), 1);

for i = 1:length(setfiles0)

    pt = setfiles0(i).folder(106:112);  % Extract subject ID
    if ismember(pt, filtered_subject_ids)

        inputfile = fullfile(setfiles0(i).folder, setfiles0(i).name);

        % Open EEGLAB
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

        % Load EEG file
        EEG = pop_loadset(inputfile);

        % Extract channel labels
        channel_labels = cell(length(EEG.chanlocs), 1);

        for c = 1:length(EEG.chanlocs)
            channel_labels{c} = EEG.chanlocs(c).labels;
        end

        % Remove channel labels that start with '_'
        channel_labels(startsWith(channel_labels, '_')) = [];

        % Store subject ID and channel labels
        subject_ids{i} = pt;
        all_channels{i} = channel_labels;
    end

end

non_empty_indices = find(~cellfun('isempty', subject_ids));

subject_ids_adol = subject_ids(~cellfun('isempty', subject_ids));
all_channels_adol = all_channels(~cellfun('isempty', all_channels));

subACWrank_adol = subACWrank(non_empty_indices);


plotInflSurf('mni', subject_ids_adol, all_channels_adol, subACWrank_adol);

saveas(gcf, [photoPath, 'allChannels_ACWColor_AdolOnly_right.png']);



%% create table with patient IDs and channel labels

% find channel rois

subject_ids = cell(length(setfiles0), 1);
all_channels = cell(length(setfiles0), 1);

for i = 1:length(setfiles0)

    pt = setfiles0(i).folder(106:112);  % Extract subject ID

    inputfile = fullfile(setfiles0(i).folder, setfiles0(i).name);

    % Open EEGLAB
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    % Load EEG file
    EEG = pop_loadset(inputfile);

    % Extract channel labels
    channel_labels = cell(length(EEG.chanlocs), 1);

    for c = 1:length(EEG.chanlocs)
        channel_labels{c} = EEG.chanlocs(c).labels;
    end

    % Remove channel labels that start with '_'
    %channel_labels(startsWith(channel_labels, '_')) = [];

    % Store subject ID and channel labels
    subject_ids{i} = pt;
    all_channels{i} = channel_labels;

end

expanded_subjects = cellfun(@(id, ch) repmat({id}, length(ch), 1), subject_ids, all_channels, 'UniformOutput', false);
expanded_subjects = vertcat(expanded_subjects{:});  % Convert to single column

% Convert channel labels to a single column
channel_labels = vertcat(all_channels{:});

% Generate channel numbers for each patient
channel_numbers = cellfun(@(ch) (1:length(ch))', all_channels, 'UniformOutput', false);
channel_numbers = vertcat(channel_numbers{:});  % Convert to single column

% Create table
T = table(expanded_subjects, channel_labels, channel_numbers, ...
          'VariableNames', {'Subject', 'Channel_Label', 'channel'});

writetable(T, '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/acw/patient_channels.csv');
