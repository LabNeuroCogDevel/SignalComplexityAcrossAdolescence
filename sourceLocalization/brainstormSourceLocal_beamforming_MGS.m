addpath resources/brainstorm3/
brainstorm nogui

%% Delay settings
delayLengths = [6 8 10];
covWindows = containers.Map([6 8 10],[5980 7980 9980]);

%% Path to FreeSurfer subjects
FSfolder = '/Volumes/Hera/preproc/7TBrainMech_rest/FS7.4.1_long'; 
FSdirs = dir(fullfile(FSfolder, '*'));

DBfolder = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/brainstorm_db';

glasserFrontalRegions = readtable('/Volumes/Hera/Projects/corticalmyelin_development/Maps/HCPMMP_glasseratlas/glasser360_regionlist_frontallobe.csv','Delimiter', '\n', 'ReadVariableNames', true);
parcelNames = extractBefore(string(glasserFrontalRegions.orig_parcelname_label), ',');

merge7t = readmatrix('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv'); 

failedSubjects = {};

%% LOOP THROUGH DELAY LENGTHS
for d = 1:length(delayLengths)

delay = delayLengths(d);
covWindow = covWindows(delay);

EEGfolder = sprintf('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/delay%d_files/', delay);

fprintf('\n===== RUNNING DELAY %d =====\n',delay)

for i = 1:length(FSdirs)

    if FSdirs(i).isdir && ~ismember(FSdirs(i).name, {'.','..'})
       
        FSsubfolder = FSdirs(i).name;

        if contains(FSsubfolder,'long')
        
            SubjectNameMRI = sprintf('%s_%s', FSsubfolder(5:9), FSsubfolder(15:22));
            subjectDBfolder = [DBfolder '/sourceLocal_' SubjectNameMRI];

            outfile = sprintf('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/glasserTimeSeries_beamforming/%s_glasserTS_delay%d.mat', SubjectNameMRI, delay);

            if ~isfile(outfile)

                try

                    sFiles = [];
                    glasserTimeSeries = {};

                    %% Protocol
                    ProtocolName = sprintf('sourceLocal_beamforming_%s', SubjectNameMRI);
                    iProtocol = bst_get('Protocol', ProtocolName);

                    if isempty(iProtocol)

                        gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
                        db_add_subject(SubjectNameMRI, [], 1);

                    else

                        gui_brainstorm('SetCurrentProtocol', iProtocol);

                    end


                    %% Delete existing delay studies
                    sStudies = bst_get('ProtocolStudies');

                    for iStudy = 1:length(sStudies.Study)

                        studyName = sStudies.Study(iStudy).Name;

                        if contains(studyName, sprintf('delay%d_files',delay))

                            fprintf('Deleting study: %s\n', studyName);
                            db_delete_studies(iStudy);

                        end

                    end


                    %% Get EEG file
                    luna = str2double(SubjectNameMRI(1:5));
                    MRIdate = str2double(SubjectNameMRI(7:14));
                    
                    rowIdx = merge7t(:,1) == luna & merge7t(:,4) == MRIdate;

                    if any(rowIdx)

                        eeg_date = merge7t(rowIdx,5);

                        pattern = sprintf('delay%d_files_%s_%s.set', delay, string(luna), string(eeg_date));
                        files = dir(fullfile(EEGfolder, pattern));

                        EEGpath = fullfile(EEGfolder, files.name);

                    else

                        continue

                    end


                    %% Import raw
                    sFiles = bst_process('CallProcess','process_import_data_epoch',[],[], ...
                        'subjectname',SubjectNameMRI, ...
                        'datafile',{EEGpath,'EEG-EEGLAB'}, ...
                        'channelreplace',1, ...
                        'channelalign',1, ...
                        'evtmode','value');


                    %% Refine registration
                    sFiles = bst_process('CallProcess','process_headpoints_refine',sFiles,[], ...
                        'tolerance',2);

                    sFiles = bst_process('CallProcess','process_channel_project',sFiles,[], ...
                        'sensortypes','EEG');


                    % BEM
                    sFiles = bst_process('CallProcess','process_generate_bem',sFiles,[], ...
                        'subjectname',SubjectNameMRI, ...
                        'nscalp',1082, ...
                        'nouter',1082, ...
                        'ninner',1082, ...
                        'thickness',4, ...
                        'method','brainstorm', ...
                        'source_abs',-1);


                    % Head model
                    sHeadModel = bst_process('CallProcess','process_headmodel',sFiles,[], ...
                        'Comment','', ...
                        'sourcespace',1, ...
                        'meg',3, ...
                        'eeg',3, ...
                        'ecog',2, ...
                        'seeg',2, ...
                        'nirs',1, ...
                        'openmeeg',struct(...
                        'BemSelect',[1,1,1], ...
                        'BemCond',[1,0.0125,1], ...
                        'BemNames',{{'Scalp','Skull','Brain'}}, ...
                        'BemFiles',{{}}, ...
                        'isAdjoint',0, ...
                        'isAdaptative',1, ...
                        'isSplit',0, ...
                        'SplitLength',4000), ...
                        'channelfile',sFiles(1).ChannelFile);


         
                    %% Covariance
                    % Process: Compute covariance (noise or data)
                    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
                        'baseline',       [], ...
                        'datatimewindow', [0, covWindow], ...
                        'sensortypes',    'MEG, EEG, SEEG, ECOG', ...
                        'target',         1, ...  % Noise covariance     (covariance over baseline time window)
                        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
                        'identity',       0, ...
                        'copycond',       0, ...
                        'copysubj',       0, ...
                        'copymatch',      0, ...
                        'replacefile',    1);  % Replace
                    
                    % Process: Compute covariance (noise or data)
                    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
                        'baseline',       [], ...
                        'datatimewindow', [0, covWindow], ...
                        'sensortypes',    'MEG, EEG, SEEG, ECOG', ...
                        'target',         2, ...  % Data covariance      (covariance over data time window)
                        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
                        'identity',       0, ...
                        'copycond',       0, ...
                        'copysubj',       0, ...
                        'copymatch',      0, ...
                        'replacefile',    1);  % Replace
                    
                    %% Sources
      
                   sSources = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
                        'output',  1, ...  % Kernel only: shared
                        'inverse', struct(...
                        'Comment',        'PNAI: EEG', ...
                        'InverseMethod',  'lcmv', ...
                        'InverseMeasure', 'nai', ...
                        'SourceOrient',   {{'free'}}, ...
                        'Loose',          0.2, ...
                        'UseDepth',       1, ...
                        'WeightExp',      0.5, ...
                        'WeightLimit',    10, ...
                        'NoiseMethod',    'median', ...
                        'NoiseReg',       0.1, ...
                        'SnrMethod',      'rms', ...
                        'SnrRms',         1e-06, ...
                        'SnrFixed',       3, ...
                        'ComputeKernel',  1, ...
                        'DataTypes',      {{'EEG'}})); 
                    
                    %% Get Glasser atlas
                    sSubject = bst_get('Subject',SubjectNameMRI);
                    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
                    CortexMat = in_tess_bst(CortexFile);
                    
                    idx = strcmp({CortexMat.Atlas.Name},'glasser');
                    glasser_atlas = CortexMat.Atlas(idx);
                    
                    ROIs = glasser_atlas.Scouts;
                    
                    
                    %% Extract scout time series
                    sScoutTS = bst_process('CallProcess','process_extract_scout',sSources,[], ...
                        'timewindow',[], ...
                        'scouts',{'glasser',{ROIs.Label}}, ...
                        'scoutfunc',1, ...
                        'isflip',0, ...
                        'isnorm',0, ...
                        'concatenate',0, ...
                        'save',1);
                    
                    
                    %% Load scout matrices
                    scoutPattern = sprintf('delay%d_*',delay);
                    
                    scoutFiles = dir(fullfile(subjectDBfolder,'data',SubjectNameMRI,scoutPattern,'matrix_scout_*.mat'));
                    
                    for s = 1:length(scoutFiles)
                        
                        glasserStruct = load([scoutFiles(s).folder '/' scoutFiles(s).name]);
                        glasserTimeSeries{s} = glasserStruct.Value;
                        
                    end
                    
                    
                    %% Save output
                    save(outfile,'glasserTimeSeries')

                catch

                    fprintf('Sub wont run: %s delay %d\n',SubjectNameMRI,delay);
                    failedSubjects{end+1} = SubjectNameMRI;

                end

            else

                fprintf('Already run: %s delay %d\n',SubjectNameMRI,delay);

            end

        else

            fprintf('Skipping folder: %s\n',FSsubfolder);

        end

    end

end

end