
%% setup env
addpath(genpath(('/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/resources/eeglab2022.1')));
addpath('/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/resources/fieldtrip-20220104');
addpath(genpath('/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/resources/EntropyHub_v2.0.0'));
ft_defaults

%% rest run one channel
inputfile = 'data/ICAwholeClean_homogenize/10129_20180919_mgs_Rem_rerefwhole_ICA_icapru.set';
EEG = read_mgs_eeg(inputfile);
delayEEG = select_delay(EEG);
channel_idx = 1;
channel_ts = delayEEG.data(channel_idx,:);
[MSx, Ci] = MSEn_channel(channel_ts);
whos
