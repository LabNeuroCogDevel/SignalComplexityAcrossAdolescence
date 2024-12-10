function ncores = ncores_from_input(inputEEG)
  [~,hostname]=system('hostname');
  
  % delay length 6 requires ~45Gb per channel
  % fix is lower and shorter
  if strncmp('rhea', hostname, 4)
     % don't care about CPU hour costs on local 'rhea'
     % and can toial away in the backround with enough memory for 10 cores.
     ncores=10;
  else
     % PSC
     % EM partition is increments of 24 cores. each increment has +1TB memory
     % ~880Gb need for 20 at once for delay len = 6
     in_type = inputEEG.event(1).type; % inputEEG.setname
     if regexp(in_type,'2')
  	   ncores=24;
     elseif  regexp(in_type,'4_6')
  	   ncores=21;
     elseif  regexp(in_type,'4_8')
  	   ncores=18;
     elseif  regexp(in_type,'4_10')
  	   ncores=16;
     else
  	   ncores=24;
     end
  end
end

%! ncores_from_input(pop_loadset('data/ICAwholeClean_homogenize/10129_20180919_mgs_Rem_rerefwhole_ICA_icapru.set'))
