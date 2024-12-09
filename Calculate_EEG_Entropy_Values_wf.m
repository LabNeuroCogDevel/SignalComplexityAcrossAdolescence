function [subjectTable] = Calculate_EEG_Entropy_Values(inputEEG, varargin)
%% multiscale entropy

% fix or delay 6, 8, 10 use different number of cores (PSC EM needed)
% always 10 on rhea
ncores = ncores_from_input(inputEEG);

% 20241118 - any errors (syntax or function not in path) within parfor will give the same cryptic error:
% otherwise get error:
%   The source code (/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/Calculate_EEG_Entropy_Values.m)
%   for the parfor-loop that is trying to execute on the worker could not be found.

t = datetime('now','TimeZone','local','Format','YYYY-MM-DD HH:mm:ss Z')
fprintf('# %s running with %d cores\n', t, ncores)
if ncores>1
  parpool('local', ncores);
  parfor c = 1:size(inputEEG.data, 1)
      t = datetime('now','TimeZone','local','Format','YYYY-MM-DD HH:mm:ss Z')
      fprintf('# %s openning per-channel cluster c=%d\n', t, c)
  
      Mobj = MSobject("SampEn");
      [MSx(c,:), Ci(:,c)] = MSEn(inputEEG.data(c,:), Mobj, 'Scales', 20, 'Methodx', 'coarse', 'RadNew', 0, 'Plotx', false);
  end
else
  for c = 1:size(inputEEG.data, 1)
      t = datetime('now','TimeZone','local','Format','YYYY-MM-DD HH:mm:ss Z')
      fprintf('# %s openning per-channel cluster c=%d\n', t, c)

      Mobj = MSobject("SampEn");
      [MSx(c,:), Ci(:,c)] = MSEn(inputEEG.data(c,:), Mobj, 'Scales', 20, 'Methodx', 'coarse', 'RadNew', 0, 'Plotx', false);
  end
end

MSxTable = array2table(MSx);
CiTable = array2table(Ci');

subjectTable = horzcat(MSxTable, CiTable);


clear Mobj
clear MSx
clear Ci


delete(gcp);
