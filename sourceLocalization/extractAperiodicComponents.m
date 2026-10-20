

%% set initial values
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/resources/brainstorm_db');

outpath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/fooof/');

%% load in all the data files
setfiles = {};
setfiles0 = dir([datapath,'/*beamforming*/data/*/1*/*specparam*.mat']);
setfiles0 = setfiles0(~startsWith({setfiles0.name}, '._'));

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(setfiles0(epo).folder, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end


subjects = cell(length(setfiles),1);

for i = 1:length(setfiles)
    subjects{i} = regexp(setfiles0(i).folder, ...
        '\d{5}_\d{8}', 'match', 'once');
end

uniqueSubjects = unique(subjects);

for s = 1:length(uniqueSubjects)

    subject = uniqueSubjects{s};

    % Files belonging to this subject
    subjectIdx = strcmp(subjects, subject);

    T = table();

    for i = find(subjectIdx)'

        subjectFile = load(setfiles{i});

        subjectExponents = subjectFile.Options.FOOOF.aperiodics;

        subTable = struct2table(subjectExponents);

        % Add subject column
        subTable.Subject = repmat({subject}, height(subTable), 1);


        % Append
        T = [T; subTable];

    end

    savePath = fullfile(outpath, ...
        [subject '_subjectAperiodicValues.csv']);

    writetable(T, savePath)

    fprintf('Saved %s (%d rows)\n', subject, height(T));

end

