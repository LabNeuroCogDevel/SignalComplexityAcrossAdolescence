
% addpath /Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/usefulCode
% addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/generalCode'))

addpath(genpath('/Users/shanemckeon/Documents/GitHub/generalCode'))
addpath(genpath('/Users/shanemckeon/Documents/GitHub/sEEG/usefulCode'))

addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/resources/eeglab2024.0'));
eeglab

addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/resources/fieldtrip-20240110'))
ft_defaults

datapath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/';
outpath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/';

resultPath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/acw_project/acw/individual_subject_files/';

% resultPath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/sEEG/acw/individual_subject_files/';

setfiles0 = dir([datapath, 'P*/rest/rest_eeglabFormat_referenced.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for i = 1:length(setfiles0)
    subject = (setfiles0(i).folder(106:112));
    clear ACWoutTable
    clear ACWout

    if ~isfile([resultPath subject '_ACW.csv'])

        inputfile = (fullfile([setfiles0(i).folder,'/', setfiles0(i).name]));

        % open eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

        EEG = pop_loadset(inputfile); % load in eeg file

        subTable = table();
        secondTrialTable = table();
        trialNumTable = table();
        channelNumTable = table();

        for t = 1:size(EEG.data, 3)

            allSecsTable = table();

            for startTime = 0:2:((EEG.times(end)/1000)-10)
                endTime = startTime+2;

                secondTrialData = pop_select(EEG, 'time', [startTime endTime], 'trial', t);

                for c = 1:size(EEG.data,1)
                    [ACWout(c,1), ACWout(c,2), acf, lags] = acw(secondTrialData.data(c,:),secondTrialData.srate, 0);

                end


                ACWoutTable = array2table(ACWout);

                secondsEpoch = repmat({sprintf('%d-%d', startTime, endTime)}, size(EEG.data, 1), 1);
                secondsEpochTable = table(secondsEpoch, 'VariableNames', {'secondsEpoch'});

                trialNum = repmat({sprintf('%d', t)}, size(EEG.data, 1), 1);
                trialNumTable = table(trialNum, 'VariableNames', {'trial'});

                channelNum = (1:size(EEG.data,1))';
                channelNumTable = table(channelNum, 'VariableNames', {'channel'});

                secondTrialTable = horzcat(ACWoutTable, channelNumTable, secondsEpochTable, trialNumTable);

                allSecsTable = vertcat(allSecsTable, secondTrialTable);

            end

            subTable = vertcat(subTable, allSecsTable);


        end

        % Create a new column with subject ID repeated for every row
        subjectIDColumn = repmat(subject, size(subTable, 1),1);

        % Add the new column to the existing table
        subTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subTable];

        subTable.Properties.VariableNames{2} = 'ACW_0';
        subTable.Properties.VariableNames{3} = 'ACW_50';

        subjectSavePath = [resultPath subject '_ACW.csv'];
        writetable(subTable, subjectSavePath)


    end
end