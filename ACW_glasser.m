
addpath((('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/code')));

%% set initial paths
datapath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/glasserTimeSeries_beamforming/');

outpath = hera('/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/individual_subject_files_beamforming/');

setfiles0 = dir([datapath,'/*glasserTS_delay10.mat']);

% remove macOS metadata files
setfiles0 = setfiles0(~startsWith({setfiles0.name}, '._'));

setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name);
end


%% Loop through all files in setfiles
for x = 1:length(setfiles)


    tokens = regexp(setfiles{x}, '(\d{5})_(\d{8})', 'tokens');

    SubjectName = tokens{1}{1};
    visitDate = tokens{1}{2};
    SubjectID = strcat(tokens{1}{1}, '_', tokens{1}{2});


    tokens = regexp(setfiles{x}, '(fix|delay8|delay6|delay10|rest)(?=\.mat$)', 'tokens');

    condition = tokens{1}{1};

    %Ensuring file does not already exist in folder

    if ~exist(['/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/individual_subject_files_beamforming/' SubjectName '_' condition 'variability.csv'], 'file')

        filename= setfiles{x};
        glasserTSstruct = load(filename);
        glasserTS = glasserTSstruct.glasserTimeSeries;

        %finding ACW using function
        for t = 1:size(glasserTS,2)
            trial = cell2mat(glasserTS(:,t));
            for c = 1:size(trial,1)
                [ACWout(c,1), ACWout(c,2), acf, lags] = acw(trial(c,:),150, 0); % Output: %   acw_0: ACW-0 %   acw_50: ACW-50
            end
            trialACW{t} = ACWout; 
        end

        % avgData = mean(cat(3, trialACW{:}), 3);
        
        trials = (cat(3, trialACW{:}));
        avgData = std(trials, 0, 3);

        ACW_0 = avgData(:,1);

        ACW_50 = avgData(:,2);

        SubjectIDTable = repmat({SubjectID}, 360, 1);

        T = table(ACW_0, ACW_50, SubjectIDTable );

        Savefile = sprintf('%s_%s_variability.csv', SubjectID, condition);

        writetable(T, fullfile(outpath, Savefile));

    end

end


% numFiles = sum(~[dir(outpath).isdir]);
% disp(numFiles);