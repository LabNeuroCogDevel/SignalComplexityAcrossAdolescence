# Preprocessing all resting state data for the Entropy project
# Pull in raw data from hera('Raw/EEG/7TBrainMech')
# run_preprocessing_pipeline.m (/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions)
# Set task as 'Resting_State'
#  Fully preprocessed data will be in Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize

matlab -nodesktop -r "addpath(genpath('../Preprocessing_Functions/')); run_preprocessing_pipeline('Resting_State')"
