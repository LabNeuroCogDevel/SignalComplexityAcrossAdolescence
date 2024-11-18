function delayEEG = select_delay(EEG)
  delayEEG = pop_selectevent( EEG, 'type',{'4_6'},'deleteevents','on');
  % Find indices of boundary events
  boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));
  
  % Remove boundary events from EEG.event structure
  delayEEG.event(boundary_indices) = [];
  
  delayEEG = pop_rmdat(delayEEG, {'4_6'},[0 5.98] ,0);
  
  % Find indices of boundary events
  boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));
  
  % Remove boundary events from EEG.event structure
  delayEEG.event(boundary_indices) = [];
end
