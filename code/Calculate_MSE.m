function Calculate_MSE(inputEEG, subject, task, resultPath, resultFolder, epoch, lengthValue, segmentLength)

if ~strcmp(epoch, 'eyesClosed') && ~strcmp(epoch, 'eyesOpen')
    if segmentLength == 0
        
        subTable = table();
        
        for t = 1:size(inputEEG.data,3)
            for c = 1:size(inputEEG.data, 1)
                Mobj = MSobject("SampEn");
                [MSx(c,:), Ci(:,c)] = MSEn(inputEEG.data(c,:,t), Mobj, 'Scales', 20, 'Methodx', 'modified', 'RadNew', 0, 'Plotx', false);
            end
            
            MSxTable = array2table(MSx);
            CiTable = array2table(Ci');
            
            trialNum = repmat({sprintf('%d', t)}, c, 1);
            trialNumTable = table(trialNum, 'VariableNames', {'trial'});
            
            channelNum = (1:64)';
            channelNumTable = table(channelNum, 'VariableNames', {'channel'});
            
            secondTrialTable = horzcat(MSxTable, CiTable,  channelNumTable, trialNumTable);
            subTable = vertcat(subTable, secondTrialTable);
            
        end
        
    elseif segmentLength == 1
        
        subTable = table();
        allSecsTable = table();
        
        for t = 1:size(inputEEG.data,3)
            % select the one second sections of the trial
            for startTime = 0:1:(lengthValue-1)
                endTime = startTime+1;
                secondTrialData = pop_select(inputEEG, 'time', [startTime endTime], 'trial', t);
                
                % run the entropy on the 1 sec that that trial
                for c = 1:size(inputEEG.data, 1)
                    Mobj = MSobject("SampEn");
                    [MSx(c,:), Ci(:,c)] = MSEn(secondTrialData.data(c,:), Mobj, 'Scales', 20, 'Methodx', 'modified', 'RadNew', 0, 'Plotx', false);
                end
                
                MSxTable = array2table(MSx);
                CiTable = array2table(Ci');
                
                secondsEpoch = repmat({sprintf('%d-%d', startTime, endTime)}, c, 1);
                secondsEpochTable = table(secondsEpoch, 'VariableNames', {'secondsEpoch'});
                
                trialNum = repmat({sprintf('%d', t)}, c, 1);
                trialNumTable = table(trialNum, 'VariableNames', {'trial'});
                
                channelNum = (1:64)';
                channelNumTable = table(channelNum, 'VariableNames', {'channel'});
                
                secondTrialTable = horzcat(MSxTable, CiTable,  channelNumTable, secondsEpochTable, trialNumTable);
                allSecsTable = vertcat(allSecsTable, secondTrialTable);
                
            end
            subTable = vertcat(subTable, allSecsTable);
        end
        
        
    elseif segmentLength == 2
        
        subTable = table();
        allSecsTable = table();
        
        for t = 1:size(inputEEG.data,3)
            % select the one second sections of the trial
            for startTime = 0:2:(lengthValue-2)
                endTime = startTime+2;
                secondTrialData = pop_select(inputEEG, 'time', [startTime endTime], 'trial', t);
                
                % run the entropy on the 1 sec that that trial
                parfor c = 1:size(inputEEG.data, 1)
                    Mobj = MSobject("SampEn");
                    [MSx(c,:), Ci(:,c)] = MSEn(secondTrialData.data(c,:), Mobj, 'Scales', 20, 'Methodx', 'modified', 'RadNew', 0, 'Plotx', false);
                end
                
                MSxTable = array2table(MSx);
                CiTable = array2table(Ci');
                
                secondsEpoch = repmat({sprintf('%d-%d', startTime, endTime)}, size(inputEEG.data, 1), 1);
                secondsEpochTable = table(secondsEpoch, 'VariableNames', {'secondsEpoch'});
                
                trialNum = repmat({sprintf('%d', t)}, size(inputEEG.data, 1), 1);
                trialNumTable = table(trialNum, 'VariableNames', {'trial'});
                
                channelNum = (1:64)';
                channelNumTable = table(channelNum, 'VariableNames', {'channel'});
                
                secondTrialTable = horzcat(MSxTable, CiTable,  channelNumTable, secondsEpochTable, trialNumTable);
                allSecsTable = vertcat(allSecsTable, secondTrialTable);
                
            end
            subTable = vertcat(subTable, allSecsTable);
        end
        
    end
    
else
    
    if segmentLength == 0
        
        subTable = table();
        
        for t = 1:size(inputEEG.data,3)
            for c = 1:size(inputEEG.data, 1)
                Mobj = MSobject("SampEn");
                [MSx(c,:), Ci(:,c)] = MSEn(inputEEG.data(c,:,t), Mobj, 'Scales', 20, 'Methodx', 'modified', 'RadNew', 0, 'Plotx', false);
            end
            
            MSxTable = array2table(MSx);
            CiTable = array2table(Ci');
            
            trialNum = repmat({sprintf('%d', t)}, c, 1);
            trialNumTable = table(trialNum, 'VariableNames', {'trial'});
            
            channelNum = (1:64)';
            channelNumTable = table(channelNum, 'VariableNames', {'channel'});
            
            secondTrialTable = horzcat(MSxTable, CiTable,  channelNumTable, trialNumTable);
            subTable = vertcat(subTable, secondTrialTable);
            
        end
        
    elseif segmentLength ==2
        
        
        subTable = table();
        allSecsTable = table();
        
        for t = 1:size(inputEEG.data,3)
            % select the one second sections of the trial
            for startTime = 0:2:((inputEEG.times(end)/1000)-10)
                endTime = startTime+2;
                secondTrialData = pop_select(inputEEG, 'time', [startTime endTime], 'trial', t);
                
                % run the entropy on the 1 sec that that trial
                parfor c = 1:size(inputEEG.data, 1)
                    Mobj = MSobject("SampEn");
                    [MSx(c,:), Ci(:,c)] = MSEn(secondTrialData.data(c,:), Mobj, 'Scales', 20, 'Methodx', 'modified', 'RadNew', 0, 'Plotx', false);
                end
                
                MSxTable = array2table(MSx);
                CiTable = array2table(Ci');
                
                secondsEpoch = repmat({sprintf('%d-%d', startTime, endTime)}, size(inputEEG.data, 1), 1);
                secondsEpochTable = table(secondsEpoch, 'VariableNames', {'secondsEpoch'});
                
                trialNum = repmat({sprintf('%d', t)}, size(inputEEG.data, 1), 1);
                trialNumTable = table(trialNum, 'VariableNames', {'trial'});
                
                channelNum = (1:64)';
                channelNumTable = table(channelNum, 'VariableNames', {'channel'});
                
                secondTrialTable = horzcat(MSxTable, CiTable,  channelNumTable, secondsEpochTable, trialNumTable);
                allSecsTable = vertcat(allSecsTable, secondTrialTable);
                
            end
            subTable = vertcat(subTable, allSecsTable);
        end
    end
    
end


% Create a new column with subject ID repeated for every row
subjectIDColumn = repmat(subject, size(subTable, 1),1);

% Add the new column to the existing table
subjectTableCols = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subTable];

subjectSavePath = [resultPath resultFolder '/individual_subject_files/' subject '_MSE_' task '_' epoch num2str(lengthValue) '.csv'];
writetable(subjectTableCols, subjectSavePath)


