
function signalComplexityCalculations(task, epoch, lengthValue, segmentLength, calculation)

addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1');
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20220104')
resultPath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/';
addpath((('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/code')));

eeglab
ft_defaults

%% Set up the folder to save results in 

if strcmp(calculation, 'ACW')

    if  strcmp(task, 'MGS')
        datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/MGS/AfterWhole/ICAwholeClean_homogenize');

        if segmentLength == 0
            resultFolder = 'acw/MGS';

        elseif segmentLength == 2
            resultFolder = 'acw/MGS/twoSecondSegments';

        end

    elseif  strcmp(task, 'Resting_State')
        datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize');

        if segmentLength == 0
            resultFolder = 'acw/Rest';

        elseif segmentLength == 2
            resultFolder = 'acw/Rest/twoSecondSegments';

        end
    end

elseif strcmp(calculation, 'MSE')

    if  strcmp(task, 'MGS')
        datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/MGS/AfterWhole/ICAwholeClean_homogenize');

        if segmentLength == 0
            resultFolder = 'mse/MGS';

        elseif segmentLength == 2
            resultFolder = 'mse/MGS/twoSecondSegments';

        end

    elseif  strcmp(task, 'Resting_State')
        datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize');

        if segmentLength == 0
            resultFolder = 'mse/Rest';

        elseif segmentLength == 2
            resultFolder = 'mse/Rest/twoSecondSegments';

        end
    end

end


%% load in all the data files

setfiles0 = dir([datapath,'/*icapru*.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues{j} = (setfiles0(j).name(1:14));
end

idvalues = unique(idvalues);
numSubj = length(idvalues);

for i = 1:numSubj
    subject = idvalues{i};
    inputfile = setfiles{i};

    EEG = pop_loadset(inputfile); % load in eeg file
    EEG = pop_eegfiltnew(EEG, 59, 61, [], 1, [], 0);

    if ~isstring(EEG.event(3).type)
        for i = 1:length(EEG.event)
            EEG.event(i).type = num2str(EEG.event(i).type);
        end
    end

    [d, currentName, ext ] = fileparts(inputfile);

    if strcmp(task, 'MGS') && strcmp(epoch, 'delay')
        if ~isfile([resultPath resultFolder '/individual_subject_files/' currentName(1:14) '_' calculation '_' task '_' epoch num2str(lengthValue) '.csv'])

            if size(EEG.data,1) > 64
                EEG = pop_select( EEG,'nochannel',{'EX3' 'EX4' 'EX5' 'EX6' 'EX7' 'EX8' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'FT7' 'FT8' 'TP7' 'TP8' 'TP9' 'TP10'});
            end

            if ~isfield(EEG.event(1), 'seconds')
                for e = 1:length(EEG.event)
                    %find event latencies in seconds
                    EEG.event(e).seconds = ([EEG.event(e).latency]-1)/EEG.srate;
                end
            end

            for e = 1:length(EEG.event)
                %find event latencies in seconds
                if e == 1
                    EEG.event(e).duration = EEG.event(e).seconds;
                elseif e == length(EEG.event)
                    EEG.event(e).duration = 2;
                else
                    EEG.event(e).duration = round(EEG.event(e+1).seconds - EEG.event(e).seconds);
                end
            end

            %select delay epochs based on length
            for e = 1:length(EEG.event)
                if num2str(EEG.event(e).type) == '4'
                    EEG.event(e).type = ['4_' num2str((EEG.event(e).duration))];
                end
            end

            try

                if lengthValue == 6
                    delayEEG = pop_selectevent( EEG, 'type',{'4_6'},'deleteevents','on');
                    % Find indices of boundary events
                    boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                    % Remove boundary events from EEG.event structure
                    delayEEG.event(boundary_indices) = [];

                    delayEEG = pop_rmdat(delayEEG, {'4_6'},[0 5.98] ,0);

                    % Find indices of boundary events
                    boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                    % Remove boundary events from EEG.event structure
                    delayEEG.event(boundary_indices) = [];

                    delayEEG = pop_epoch(delayEEG, {'4_6'}, [0 5.98], 'epochinfo', 'yes'); % create epochs using selected events

                elseif lengthValue == 8
                    delayEEG = pop_selectevent( EEG, 'type',{'4_8'},'deleteevents','on');
                    % Find indices of boundary events
                    boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                    % Remove boundary events from EEG.event structure
                    delayEEG.event(boundary_indices) = [];

                    delayEEG = pop_rmdat(delayEEG, {'4_8'},[0 7.98] ,0);

                    % Find indices of boundary events
                    boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                    % Remove boundary events from EEG.event structure
                    delayEEG.event(boundary_indices) = [];

                    delayEEG = pop_epoch( delayEEG, {'4_8'}, [0 7.98], 'epochinfo', 'yes'); % create epochs using selected events

                elseif lengthValue == 10
                    delayEEG = pop_selectevent( EEG, 'type',{'4_10'},'deleteevents','on');
                    % Find indices of boundary events
                    boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                    % Remove boundary events from EEG.event structure
                    delayEEG.event(boundary_indices) = [];

                    delayEEG = pop_rmdat(delayEEG, {'4_10'},[0 9.98] ,0); % Find indices of boundary events

                    boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                    % Remove boundary events from EEG.event structure
                    delayEEG.event(boundary_indices) = [];

                    delayEEG = pop_epoch(delayEEG, {'4_10'}, [0 9.98], 'epochinfo', 'yes'); % create epochs using selected event

                end

                if (length(delayEEG.event) == length(EEG.event)) || (length(delayEEG.event) < 50 && lengthValue == 6) || (length(delayEEG.event) < 20 && lengthValue == 8)|| (length(delayEEG.event) < 5 && lengthValue == 10)
                    disp("too few triggers")
                    continue; % Move to next person in the loop if no events were removed
                end

                if strcmp(calculation, 'ACW')
                    calculate_ACW(delayEEG, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength)
                elseif strcmp(calculation, 'MSE')
                    Calculate_MSE(delayEEG, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength);
                end
                
            catch
                disp("did not run")
                continue;
            end


        end

    elseif strcmp(task, 'MGS') && strcmp(epoch, 'fix')
        if ~isfile([resultPath resultFolder '/individual_subject_files/' currentName(1:14) '_' calculation '_' task '_' epoch num2str(lengthValue) '.csv'])

            if size(EEG.data,1) > 64
                EEG = pop_select( EEG,'nochannel',{'EX3' 'EX4' 'EX5' 'EX6' 'EX7' 'EX8' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'FT7' 'FT8' 'TP7' 'TP8' 'TP9' 'TP10'});
            end


            try
                fixEEG = pop_selectevent( EEG, 'type',{'2'},'deleteevents','on');
                % Find indices of boundary events
                boundary_indices = find(strcmp({fixEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                fixEEG.event(boundary_indices) = [];

                fixEEG = pop_rmdat(fixEEG, {'2'},[0 1.98] ,0);

                % Find indices of boundary events
                boundary_indices = find(strcmp({fixEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                fixEEG.event(boundary_indices) = [];

                fixEEG = pop_epoch(fixEEG, {'2'}, [0 1.98], 'epochinfo', 'yes'); % create epochs using selected event


                if (length(fixEEG.event) == length(EEG.event)) || (length(fixEEG.event) < 80)
                    disp("too few triggers")
                    continue; % Move to next person in the loop if no events were removed
                end

                if strcmp(calculation, 'ACW')
                    calculate_ACW(fixEEG, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength)
                elseif strcmp(calculation, 'MSE')
                    Calculate_MSE(fixEEG, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength);
                end

            catch
                disp("did not run")
                continue;
            end


        end

    elseif strcmp(task, 'Resting_State')
        if strcmp(epoch, 'eyesClosed')

            if ~isfile([resultPath resultFolder '/individual_subject_files/' currentName(1:14) '_' calculation '_rest_eyesClosed.csv'])

                EEGclosedeyes = pop_rmdat(EEG, {'16129', '15261','15361','0'},[0 4] ,0);

                % Find indices of boundary events
                boundary_indices = find(strcmp({EEGclosedeyes.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                EEGclosedeyes.event(boundary_indices) = [];

                if (length(EEGclosedeyes.event) == length(EEG.event)) || (length(EEGclosedeyes.event) < 40) || (length(EEGclosedeyes.event) > 100)
                    disp("wrong nunmber of triggers")
                    continue; % Move to next person in the loop if no events were removed
                end
               

                if strcmp(calculation, 'ACW')
                    calculate_ACW(EEGclosedeyes, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength)
                elseif strcmp(calculation, 'MSE')                  
                    Calculate_MSE(EEGclosedeyes, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength);
                end
                
            end

        elseif strcmp(epoch, 'eyesOpen')

            if ~isfile([resultPath resultFolder '/individual_subject_files/' currentName(1:14) '_' calculation '_' task '_' epoch num2str(lengthValue) '.csv'])

                EEGopeneyes = pop_rmdat(EEG, {'16130', '15362','1'},[0 4] ,0);

                % Find indices of boundary events
                boundary_indices = find(strcmp({EEGopeneyes.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                EEGopeneyes.event(boundary_indices) = [];

                if (length(EEGopeneyes.event) == length(EEG.event)) || (length(EEGopeneyes.event) < 40) || (length(EEGopeneyes.event) > 100)
                    disp("wrong nunmber of triggers")
                    continue; % Move to next person in the loop if no events were removed
                end

                if strcmp(calculation, 'ACW')
                    calculate_ACW(EEGopeneyes, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength)
                elseif strcmp(calculation, 'MSE')
                    Calculate_MSE(EEGopeneyes, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength);
                end
            end

        end
    end

end
end















