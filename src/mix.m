function mixed_signals = mix(X, mixture_type='linear')
  n = size(X)(1);

  if strcmp(mixture_type, 'conv')
    'CONV'
    mixed_signals = 'CONV';
  else
    A = unifrnd(1, 5, n, n);
    S = X' * A;
    mixed_signals = S';
  end
end
