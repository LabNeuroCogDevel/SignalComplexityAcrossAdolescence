function [subjectTable] = Calculate_EEG_Entropy_Values(inputEEG, varargin)

%% multiscale entropy

[~,hostname]=system('hostname');
if strncmp('rhea',hostname,4)
   ncores=10; % 64
   parpool('local', ncores);
else
   % PSC
   ncores=20; % 64
   parpool('local', ncores);
end


% 20241118 - parpool needs to know about this file
% otherwise get error:
%   The source code (/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/Calculate_EEG_Entropy_Values.m)
%   for the parfor-loop that is trying to execute on the worker could not be found.
%addAttachedFiles(gcp,["/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/Calculate_EEG_Entropy_Values.m",...
%                      "resources/" ])

% pool = gcp();
%this_script='/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/Calculate_EEG_Entropy_Values.m';
%addAttachedFiles(pool, this_script);


parfor c = 1:size(inputEEG.data, 1)
    fprintf('# openning per-channel cluster c=%d\n',c)
    % 20241118 - parpool complaining about this file (current one) being MIA
    % disp(dir(this_script))

    Mobj = MSobject("SampEn");
    [MSx(c,:), Ci(:,c)] = MSEn(inputEEG.data(c,:), Mobj, 'Scales', 20, 'Methodx', 'coarse', 'RadNew', 0, 'Plotx', false);

end

MSxTable = array2table(MSx);
CiTable = array2table(Ci');

subjectTable = horzcat(MSxTable, CiTable);


clear Mobj
clear MSx
clear Ci


delete(gcp);

%% spectral entropy

% for j = 1:length(setfiles0)
%
%     idvalues(j,:) = (setfiles0(j).name(1:14));
%     inputfile = setfiles{j};
%     if ~isfile([entropyPath idvalues(j,:) '_SpectralEntropy_broadband.csv'])
%
%         EEG = pop_loadset(inputfile); % load in eeg file
%         EEGopeneyes = pop_rmdat(EEG, {'16130'},[0 4] ,0);
%         onemin = EEGopeneyes.data(:,1:9000);
%
%
%         parpool('local', 25);
%
%         parfor c = 1:size(onemin, 1)
%
%             [Spec(c,:), gammaBandEn(c,:)] = SpecEn(onemin(c,:), 'N', 150, 'Freqs', [.4, 1], 'Logx', exp(1), 'Norm' , true);
%             [Spec(c,:), betaBandEn(c,:)] = SpecEn(onemin(c,:), 'N', 150, 'Freqs', [.16, .4], 'Logx', exp(1), 'Norm' , true);
%             [Spec(c,:), alphaBandEn(c,:)] = SpecEn(onemin(c,:), 'N', 150, 'Freqs', [.1, .16], 'Logx', exp(1), 'Norm' , true);
%             [Spec(c,:), thetaBandEn(c,:)] = SpecEn(onemin(c,:), 'N', 150, 'Freqs', [.04, .1], 'Logx', exp(1), 'Norm' , true);
%
%         end
%         delete(gcp('nocreate'));
%
%
%         gammaBandEnTable = array2table(gammaBandEn);
%         betaBandEnTable = array2table(betaBandEn);
%         alphaBandEnTable = array2table(alphaBandEn);
%         thetaBandEnTable = array2table(thetaBandEn);
%         SpecTable = array2table(Spec);
%
%
%
%         subjectTable = horzcat(gammaBandEnTable, betaBandEnTable, alphaBandEnTable, thetaBandEnTable, SpecTable);
%
%         % Create a new column with subject ID repeated for every row
%         subjectIDColumn = repmat(idvalues(j,:), size(subjectTable, 1),1);
%
%         % Add the new column to the existing table
%         subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];
%
%
%         subjectSavePath = [entropyPath idvalues(j,:) '_SpectralEntropy_broadband.csv'];
%         writetable(subjectTable, subjectSavePath)
%
%
%     end
%
%     clear gammaBandEn
%     clear betaBandEn
%     clear alphaBandEn
%     clear thetaBandEn
%     clear Spec
%
% end
