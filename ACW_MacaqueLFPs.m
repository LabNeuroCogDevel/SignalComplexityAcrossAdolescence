

function playingWithLFPs(task)

addpath((('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/')));
addpath((('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/code')));
addpath('/Volumes/Hera/Shane/LFPData-20250428T155619Z-001/LFPData/')


if strcmp(task, 'delay')

    resultPath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MonkeyData/individual_files_delay';

elseif strcmp(task, 'fix')

    resultPath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/MonkeyData/individual_files_fix';

end


channel_info = readtable('channel_info.xlsx');

get_filename = @(row) sprintf('%s_CH%d.mat', row.session{1}, row.channel); % Naming convention made up of behavior session name and electrode number
% class_to_radian = @(class) 2 * pi * (class - 1)/8; % classses increase in radian 1 through 8 starting at 0.
cue_duration = 0.5; % seconds

for i_electrode=1:size(channel_info,1)
    
    delay_duration = channel_info.delay_duration(i_electrode); % Time between cue offset and go signal
    load(get_filename(channel_info(i_electrode, :)), 'LFPData');
    
    filename = get_filename(channel_info(i_electrode, :));
    monkey = filename(1:3);
    behavSess = filename(4:6);
    channel = filename(10:12);
    
    if ~isfile([resultPath '/' monkey '_acw_' 'behavSess' behavSess '_' channel '.csv'])
        
        
        get_lfp_slice = @(fs, t0, t_slice) [round(t_slice(1)*fs):round(t_slice(2)*fs)] + round(t0*fs);
   
        fs = LFPData.class(1).ntr(1).Sampling_Rate;
         
        if strcmp(task, 'delay')

            target_time = [cue_duration, cue_duration+delay_duration]; % Time around cue onset from which LFP is to be aligned to

        elseif strcmp(task, 'fix')
            
            target_time = [-1, 0]; % Time around cue onset from which LFP is to be aligned to

        end

        dummy_lfp_slice = get_lfp_slice(fs, 0, target_time);
        t_slice = dummy_lfp_slice/fs;
        
        clear lfp_current_class
        
        lfp_by_class = cell(numel(LFPData.class), 1);  % one cell per class
        
        for i_class = 1:numel(LFPData.class)
            n_trials = numel(LFPData.class(i_class).ntr);
            lfp_trials = zeros(n_trials, numel(get_lfp_slice(fs, 0, target_time)));  % preallocate per class
            
            for i_trial = 1:n_trials
                t_idx = get_lfp_slice(fs, LFPData.class(i_class).ntr(i_trial).Cue_onT, target_time);
                try
                lfp_trials(i_trial, :) = LFPData.class(i_class).ntr(i_trial).LFP(t_idx);
                catch
                end
            end
            
            lfp_by_class{i_class} = lfp_trials;  % store trials (2D array) for each class
        end
        
        lfp_all_trials = vertcat(lfp_by_class{:});  % size: [total_trials x timepoints]

        
        % plot(t_slice, mean(lfp_all_trials,1))
        % plot(t_slice, lfp_all_trials(2,:))
        
        clear ACWout
        
        for t = 1:size(lfp_all_trials,1)
            lfp_down = resample(lfp_all_trials(t,:), 150, 500);

            ACWout{t,1} = monkey;
            ACWout{t,2} = behavSess;
            ACWout{t,3} = channel;
            ACWout{t,4} = t;
            [ACWout{t,5}, ACWout{t,6}, acf, lags] = acw(lfp_down,150, 1);
        end
        
        % Add the new column to the existing table
        subjectTable = array2table(ACWout, 'VariableNames', {'Monkey', 'BehavSess', 'Channel', 'Trial','ACW_0', 'ACW_50'});
        
        subjectTable.NeuronArea = repmat(cell2mat(channel_info(i_electrode, :).Neuron_area), size(ACWout,1),1);
        subjectTable.MatureAge = repmat((channel_info(i_electrode, :).mature_age), size(ACWout,1),1);
        subjectTable.chronological_age = repmat((channel_info(i_electrode, :).Neuron_age + channel_info(i_electrode, :).mature_age), size(ACWout,1),1);
        subjectTable.group = repmat((channel_info(i_electrode, :).group), size(ACWout,1),1);
        subjectTable.mature_group = repmat((channel_info(i_electrode, :).mature_group), size(ACWout,1),1);
        subjectTable.age_group = repmat((channel_info(i_electrode, :).age_group), size(ACWout,1),1);

        subjectSavePath = [resultPath '/' monkey '_acw_' 'behavSess' behavSess '_' channel '.csv'];
        writetable(subjectTable, subjectSavePath)
        
    end
    
end

