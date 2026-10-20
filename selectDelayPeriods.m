function selectDelayPeriods()

addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1');
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20220104')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/code');

eeglab
ft_defaults

%% Paths
datapath  = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/MGS/AfterWhole/ICAwholeClean_homogenize';
resultpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/';

%% Get files
setfiles0 = dir(fullfile(datapath,'*icapru*.set'));
setfiles  = fullfile(datapath,{setfiles0.name});

idvalues = unique(cellfun(@(x) x(1:14), {setfiles0.name}, 'UniformOutput', false));
numSubj  = length(idvalues);

%% Delay lengths to process
delayLengths = [6 8];

for i = 1:numSubj
    
    inputfile = setfiles{i};
    [~, currentName, ~] = fileparts(inputfile);
    nameShort = currentName(1:min(14,length(currentName)));

    try
        EEG = pop_loadset(inputfile);
        EEG = pop_eegfiltnew(EEG, 59, 61, [], 1, [], 0);

        %% Ensure event types are strings
        for e = 1:length(EEG.event)
            if ~ischar(EEG.event(e).type)
                EEG.event(e).type = num2str(EEG.event(e).type);
            end
        end

        %% Compute seconds + duration
        for e = 1:length(EEG.event)
            EEG.event(e).seconds = (EEG.event(e).latency-1)/EEG.srate;
        end

        for e = 1:length(EEG.event)
            if e == 1
                EEG.event(e).duration = EEG.event(e).seconds;
            elseif e == length(EEG.event)
                EEG.event(e).duration = 2;
            else
                EEG.event(e).duration = round(EEG.event(e+1).seconds - EEG.event(e).seconds);
            end
        end

        %% Label delay events
        for e = 1:length(EEG.event)
            if strcmp(EEG.event(e).type,'4')
                EEG.event(e).type = ['4_' num2str(EEG.event(e).duration)];
            end
        end

        %% ---- LOOP THROUGH 6, 8, 10 ----
        for L = delayLengths
            
            resultFolder = ['delay' num2str(L) '_files'];
            savePath = fullfile(resultpath,resultFolder);

            if ~exist(savePath,'dir')
                mkdir(savePath);
            end

            saveName = [resultFolder '_' nameShort '.set'];
            
            if ~isfile(fullfile(savePath,saveName))

                epochEnd = L - 0.02;   % 5.98, 7.98, 9.98
                eventType = ['4_' num2str(L)];

                delayEEG = pop_selectevent(EEG,'type',{eventType},'deleteevents','on');

                % Remove boundaries
                boundary_idx = strcmp({delayEEG.event.type},'boundary');
                delayEEG.event(boundary_idx) = [];

                delayEEG = pop_rmdat(delayEEG,{eventType},[0 epochEnd],0);

                % Remove boundaries again
                boundary_idx = strcmp({delayEEG.event.type},'boundary');
                delayEEG.event(boundary_idx) = [];

                delayEEG = pop_epoch(delayEEG,{eventType},[0 epochEnd],'epochinfo','yes');

                delayEEG = pop_saveset(delayEEG,...
                    'filename',saveName,...
                    'filepath',savePath);
            end
        end

    catch
        disp(['Did not run: ' nameShort])
        continue
    end
end

end