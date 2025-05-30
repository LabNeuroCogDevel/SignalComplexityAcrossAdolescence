
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/usefulCode'));
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/generalCode'));
datapath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/';
outpath = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/';

%load in all the files
setfiles0 = dir([datapath, 'P*/rest/rest*ns2*']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1:length(setfiles0)

    pt = (setfiles0(j).folder(106:112));

    if isfile(['/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/sEEG_backup/sEEG_rawData/', pt, '/rest/rest_eeglabFormat_referenced.set']) 
    else

        [matFile, info] = sEEG_ns_to_mat(pt, 'rest', setfiles0(j));

        setfile = sEEG_mat_to_set(pt, 'rest' , matFile, info);

    end

end