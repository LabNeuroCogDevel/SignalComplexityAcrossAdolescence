function [MSx, Ci] = MSEn_channel(channel_ts, varargin)
    % input_data = inputEEG.data(c,:);
    Mobj = MSobject("SampEn");
    [MSx, Ci] = MSEn(channel_ts, Mobj, 'Scales', 20, 'Methodx', 'coarse', 'RadNew', 0, 'Plotx', false);
end
