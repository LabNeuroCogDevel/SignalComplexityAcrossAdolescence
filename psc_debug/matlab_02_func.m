function res = matlab_02_func(data, n_parts)
  % reuse existing pool. for debugging: in prod func only called once per matlab instance/node
  p = gcp('nocreate');
  if isempty(p), p=parpool(2); end
  n = size(data,2);
  res = zeros(n, n_parts);
  fprintf('running %d datapoints through  %d iterations on %d workers\n',size(res), p.NumWorkers);
  parfor i = 1:n_parts
     fprintf('  i=%d\n',i);
     res(:,i) = rand(1,n).*data;
  end
end
