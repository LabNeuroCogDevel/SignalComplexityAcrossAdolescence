
function calculate_ACW(inputEEG, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength)
if ~strcmp(epoch, 'eyesClosed') && ~strcmp(epoch, 'eyesOpen')


    if segmentLength == 0

        acwChan=[];

        for c = 1:size(inputEEG.data,1)
            for t = 1:size(inputEEG.data, 3)
                ACWout{t,1} = inputEEG.chanlocs(c).labels;
                ACWout{t,2} = t;
                [ACWout{t,3}, ACWout{t,4}, acf{c,t}, lags{c,t}] = acw(inputEEG.data(c,:,t),inputEEG.srate, 0);
            end
            acwChan = [acwChan; ACWout];
        end

        % Create a new column with subject ID repeated for every row
        subjectIDColumn = repmat(subject, size(acwChan, 1),1);

        % Add the new column to the existing table
        subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), acwChan];
        subjectTable.Properties.VariableNames{2} = 'Channel';
        subjectTable.Properties.VariableNames{3} = 'Trial';
        subjectTable.Properties.VariableNames{4} = 'ACW_0';
        subjectTable.Properties.VariableNames{5} = 'ACW_50';

        subjectSavePath = [resultPath resultFolder '/individual_subject_files/' subject '_acw_' task '_' epoch num2str(lengthValue) '.csv'];
        writetable(subjectTable, subjectSavePath)

        subjectSavePathACF = [resultPath resultFolder '/individual_subject_files/' subject '_ACF_' task '_' epoch num2str(lengthValue) '.mat'];
        save(subjectSavePathACF, 'acf')

        subjectSavePathlags = [resultPath resultFolder '/individual_subject_files/' subject '_lags_' task '_' epoch num2str(lengthValue) '.mat'];
        save(subjectSavePathlags, 'lags')


    elseif segmentLength == 2

        subTable = table();
        results = struct();


        for t = 1:size(inputEEG.data, 3)
            allSecsTable = table();


            for startTime = 0:2:(lengthValue-2)
                endTime = startTime+2;

                secondTrialData = pop_select(inputEEG, 'time', [startTime endTime], 'trial', t);

                for c = 1:size(inputEEG.data,1)
                    [ACWout(c,1), ACWout(c,2), acf, lags] = acw(secondTrialData.data(c,:),secondTrialData.srate, 0);

                end


                ACWoutTable = array2table(ACWout);

                secondsEpoch = repmat({sprintf('%d-%d', startTime, endTime)}, size(inputEEG.data, 1), 1);
                secondsEpochTable = table(secondsEpoch, 'VariableNames', {'secondsEpoch'});

                trialNum = repmat({sprintf('%d', t)}, size(inputEEG.data, 1), 1);
                trialNumTable = table(trialNum, 'VariableNames', {'trial'});

                channelNum = (1:64)';
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

        subjectSavePath = [resultPath resultFolder '/individual_subject_files/' subject '_ACW_' task '_' epoch num2str(lengthValue) '.csv'];
        writetable(subTable, subjectSavePath)


    end

else

    if segmentLength == 2

        subTable = table();
        results = struct();


        for t = 1:size(inputEEG.data, 3)
            allSecsTable = table();


            for startTime = 0:2:((inputEEG.times(end)/1000)-10)
                endTime = startTime+2;

                secondTrialData = pop_select(inputEEG, 'time', [startTime endTime], 'trial', t);

                for c = 1:size(inputEEG.data,1)
                    [ACWout(c,1), ACWout(c,2), acf, lags] = acw(secondTrialData.data(c,:),secondTrialData.srate, 0);

                end


                ACWoutTable = array2table(ACWout);

                secondsEpoch = repmat({sprintf('%d-%d', startTime, endTime)}, size(inputEEG.data, 1), 1);
                secondsEpochTable = table(secondsEpoch, 'VariableNames', {'secondsEpoch'});

                trialNum = repmat({sprintf('%d', t)}, size(inputEEG.data, 1), 1);
                trialNumTable = table(trialNum, 'VariableNames', {'trial'});

                channelNum = (1:64)';
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

        subjectSavePath = [resultPath resultFolder '/individual_subject_files/' subject '_ACW_' task '_' epoch num2str(lengthValue) '.csv'];
        writetable(subTable, subjectSavePath)
    end




end


