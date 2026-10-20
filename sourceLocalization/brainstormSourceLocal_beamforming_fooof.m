addpath resources/brainstorm3/
brainstorm nogui


%% Path to FreeSurfer subjects
FSfolder = '/Volumes/Hera/preproc/7TBrainMech_rest/FS7.4.1_long';
FSdirs = dir(fullfile(FSfolder, '*'));

EEGfolder = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize/';
DBfolder = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/brainstorm_db';
glasserFrontalRegions = readtable('/Volumes/Hera/Projects/corticalmyelin_development/Maps/HCPMMP_glasseratlas/glasser360_regionlist_frontallobe.csv','Delimiter', '\n', 'ReadVariableNames', true);
parcelNames = extractBefore(string(glasserFrontalRegions.orig_parcelname_label), ',');
merge7t = readmatrix('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv');
%
failedSubjects = {};

for i = 1:length(FSdirs)
    if FSdirs(i).isdir && ~ismember(FSdirs(i).name, {'.', '..'})
        
        FSsubfolder = FSdirs(i).name;  % folder name
        
        
        % Only include folders with an underscore (ID_date format)
        if contains(FSsubfolder, 'long')
            
            tokens = regexp(FSsubfolder, 'sub-(\d{5})_ses-(\d{8})', 'tokens');

            SubjectNameMRI = tokens{1}{1};
            visitDate = tokens{1}{2};
            SubjectID = strcat(tokens{1}{1}, '_', tokens{1}{2});

            subjectDBfolder = [DBfolder '/sourceLocal_beamforming_' SubjectID];
            
            %% delete previous fooof - if need be
            
            
            searchPattern = fullfile(subjectDBfolder, 'data', SubjectID, '**', '*timefreq*.mat');
            files = dir(searchPattern);
            
            for i = 1:numel(files)
                fullFile = fullfile(files(i).folder, files(i).name);
                fprintf('Deleting %s\n', fullFile);
                delete(fullFile);
            end
            
            %%
  
            
            searchPattern = fullfile(subjectDBfolder, 'data', SubjectID, '**', '*specparam.mat');
            files = dir(searchPattern);
            
            if isempty(files)
  
                try
                    sFiles = [];
                    
                    % Set a unique protocol for each subject
                    ProtocolName = sprintf('sourceLocal_beamforming_%s', SubjectID);
                    
                    % Get existing protocols
                    iProtocol = bst_get('Protocol', ProtocolName);
                    if isempty(iProtocol)
                        
                        % Create a new protocol for each subject
                        gui_brainstorm('CreateProtocol', ProtocolName, 0, 0); %
                        
                        % Add the subject to the protocol (if not already present)
                        db_add_subject(SubjectNameMRI, []);
                    else
                        gui_brainstorm('SetCurrentProtocol', iProtocol);
                        
                    end
                    
                    
                  
                    scoutFiles = dir(fullfile(subjectDBfolder,'data',SubjectID,'1*','matrix_scout_*.mat'));

                    sFiles = cell(length(scoutFiles),1);

                    for i = 1:length(scoutFiles)
                        fullpath = fullfile(scoutFiles(i).folder, scoutFiles(i).name);

                        % Remove everything before the Brainstorm data directory
                        idx = strfind(fullpath, [filesep 'data' filesep]);

                        sFiles{i} = fullpath(idx + length([filesep 'data' filesep]):end);
                    end

                    bst_report('Start', sFiles);


                    % Process: Power spectrum density (Welch)
                    sPSD = bst_process('CallProcess', 'process_psd', sFiles, [], ...
                        'timewindow',  [], ...
                        'win_length',  2, ...
                        'win_overlap', 50, ...
                        'units',       'physical', ...  % Physical: U2/Hz
                        'sensortypes', 'EEG', ...
                        'win_std',     0, ...
                        'edit',        struct(...
                        'Comment',         'Power', ...
                        'TimeBands',       [], ...
                        'Freqs',           [], ...
                        'ClusterFuncTime', 'none', ...
                        'Measure',         'power', ...
                        'Output',          'all', ...
                        'SaveKernel',      0));

                    
                    % Process: specparam: Fitting oscillations and 1/f
                    sFOOOF = bst_process('CallProcess', 'process_fooof', sPSD, [], ...
                        'implementation', 'matlab', ...  % Matlab
                        'freqrange',      [1, 60], ...
                        'powerline',      '60', ...  % 60 Hz
                        'method',         'leastsquare', ...  % Default
                        'peakwidth',      [0.5, 12], ...
                        'maxpeaks',       3, ...
                        'minpeakheight',  3, ...
                        'proxthresh',     2, ...
                        'apermode',       'fixed', ...  % Fixed
                        'guessweight',    'none', ...  % None
                        'sorttype',       'param', ...  % Peak parameters
                        'sortparam',      'frequency', ...  % Frequency
                        'sortbands',      {'delta', '2, 4'; 'theta', '5, 7'; 'alpha', '8, 12'; 'beta', '15, 29'; 'gamma1', '30, 59'; 'gamma2', '60, 90'});


                catch
                    fprintf('Sub wont run: %s\n', SubjectNameMRI);
                    failedSubjects{end+1} = SubjectNameMRI;
                    
                end
                
            else
                fprintf('Sub already run, skipping %s\n', SubjectNameMRI);
                
            end
            
        else
            fprintf('Skipping folder (no date): %s\n', FSsubfolder);
            
        end
    end
    
    
end
