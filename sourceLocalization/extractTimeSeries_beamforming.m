%% Path to FreeSurfer subjects

addpath resources/brainstorm3/

brainstorm nogui

addpath /Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sourceLocalization/

DBfolder = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/brainstorm_db';
bst_set('BrainstormDbDir', DBfolder);

FSdirs = dir(fullfile(DBfolder, '*beamforming*'));
glasserRegions = readtable('/Volumes/Hera/Projects/corticalmyelin_development/Maps/HCPMMP_glasseratlas/glasser360_regionlist.csv','Delimiter', '\n', 'ReadVariableNames', true);
parcelNames = extractBefore(string(glasserRegions.orig_parcelname_label), ',');

merge7t = readmatrix('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv');

for i = 1:length(FSdirs)
    
    subfolder = FSdirs(i).name;  % folder name
    
    SubjectName = sprintf('%s_%s', subfolder(25:29), subfolder(31:38));
    
    subjectDatafolder = [subfolder '/data/' SubjectName];
    
    searchPattern = fullfile(DBfolder, subjectDatafolder, '**', 'data_concat_*.mat');
    
    files = dir(searchPattern);
    
    if ~exist(['/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/glasserTimeSeries_beamforming/' SubjectName '_glasserTS.mat'], 'file')
        
        if ~isempty(files)
            
            
            sFiles = [];
            
            % Name of the protocol
            ProtocolName = sprintf('sourceLocal_beamforming_%s', SubjectName);
            
            % Check if the protocol already exists
            iProtocol = bst_get('Protocol', ProtocolName);
            
            if isempty(iProtocol)
                % Create the protocol (only happens first time)
                gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
                
                % Add subject without GUI warnings (overwrite=1)
                db_add_subject(SubjectNameMRI, [], 1);
            else
                % Just open the protocol and continue working
                gui_brainstorm('SetCurrentProtocol', iProtocol);
            end
            
            
            sSubject = bst_get('Subject', SubjectName);
            CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
            CortexMat = in_tess_bst(CortexFile);
            
            idx = strcmp({CortexMat.Atlas.Name}, 'glasser');
            glasser_atlas = CortexMat.Atlas(idx);
            % Extract all labels from the Glasser atlas
            labels = {glasser_atlas.Scouts.Label};
            
            % Remove the trailing ' L' or ' R' for matching
            cleanLabels = regexprep(labels, ' [LR]$', '');
            
            % Find matches
            ROIsIdx = find(ismember(cleanLabels, parcelNames));
            
            ROIs = glasser_atlas.Scouts(ROIsIdx);
            
            searchPatternBeam = fullfile(DBfolder, subjectDatafolder, '**', 'results_PNAI_EEG_KERNEL_*.mat');
            
            filesBeam = dir(searchPatternBeam);
            
            %'link|10129_20180917/16130/results_sLORETA_EEG_KERNEL_251022_1405.mat|10129_20180917/16130/data_concat_251022_1352.mat'}
            %format
            
            if ~isempty(filesBeam)
                [~, lastFolder] = fileparts(filesBeam(end).folder);
                
                
                sFiles = {['link|' SubjectName '/' lastFolder '/' filesBeam(end).name '|' SubjectName '/' lastFolder '/' files(end).name]};
                
                % Start a new report
                bst_report('Start', sFiles);
                
                
                
                %extrsact scout time series
                sScoutTS = bst_process('CallProcess', 'process_extract_scout', sFiles, [], ...
                    'timewindow',   [], ...                  % full time window
                    'scouts',       {'glasser', {ROIs.Label}}, ...  % use the atlas name + scout labels
                    'scoutfunc',    1, ...                   % 1 = mean, 2 = max, 3 = PCA, etc.
                    'isflip',       0, ...                   % don't flip sign
                    'isnorm',       0, ...                   % no normalization
                    'concatenate',  1, ...                    % concatenate all scouts
                    'save',         1);                      % save results in the database
                
                searchglasser360Scouts = fullfile(DBfolder, subjectDatafolder(1:44), sScoutTS.FileName);
                
                filesGlasser = dir(searchglasser360Scouts);
                
                
                glasserTimeSeries = load([filesGlasser.folder '/' filesGlasser.name]);
                
                glasserTS = glasserTimeSeries.Value;
                save(['/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/glasser/glasserTimeSeries_beamforming/' SubjectName '_glasserTS.mat'], 'glasserTS')
                
                
                %                    glasserRegionOrder =  glasserTimeSeries.Description;
                %                    save(['/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/glasserTSorder.mat'], 'glasserRegionOrder')
            else
                continue
            end
            
            
            
            
            
            
        end
    end
end

