function matlab_01_entry(n_parts)
  data = 1:10;
  fprintf('# starting matlab_02_func\n')
  % 01_sbatch.bash: 4 parts, matlab_02_func: 2 pools
  x = matlab_02_func(data, n_parts);
  fprintf('# finished matlab_02_func\n')
  disp(sum(x));
end
