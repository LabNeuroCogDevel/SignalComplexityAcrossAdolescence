function EEG = read_mgs_eeg(inputfile)

  EEG = pop_loadset(inputfile); % load in eeg file
  [d, currentName, ext ] = fileparts(inputfile);
  if size(EEG.data,1) > 64
      EEG = pop_select( EEG,'nochannel',{'EX3' 'EX4' 'EX5' 'EX6' 'EX7' 'EX8' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'FT7' 'FT8' 'TP7' 'TP8' 'TP9' 'TP10'});
  end

  if ~isfield(EEG.event(1), 'seconds')
      for e = 1:length(EEG.event)
          %find event latencies in seconds
          EEG.event(e).seconds = ([EEG.event(e).latency]-1)/EEG.srate;
      end
  end

  for e = 1:length(EEG.event)
      %find event latencies in seconds
      if e == 1
          EEG.event(e).duration = EEG.event(e).seconds;
      elseif e == length(EEG.event)
          EEG.event(e).duration = 2;
      else
          EEG.event(e).duration = round(EEG.event(e+1).seconds - EEG.event(e).seconds);
      end
  end

  %select delay epochs based on length
  for e = 1:length(EEG.event)
      if num2str(EEG.event(e).type) == '4'
          EEG.event(e).type = ['4_' num2str((EEG.event(e).duration))];
      end
  end
end
