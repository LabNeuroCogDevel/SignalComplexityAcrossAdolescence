function res = matlab_02_func(data, n_parts)
  parpool(2)
  n = size(data,1)
  res = zeros(n, n_parts)
  parfor i = 1:n_parts
     res(:,i) = rand(n,1)*data;
  end
end
