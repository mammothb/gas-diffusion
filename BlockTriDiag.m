function T = BlockTriDiag(n, Tmain, Tsub, Tsup)
  % verify the inputs in this mode are 2-d arrays.
  if length(size(Tmain)) ~= 2 || length(size(Tsub)) ~= 2 ||...
      length(size(Tsup)) ~= 2
    error('Inputs must be 2d arrays if a replication factor is provided');
  end
  % get block size and check for consistency
  [p, q] = size(Tmain);
  if isempty(Tmain)
    error('Blocks must be non-empty arrays or scalars');
  end
  if any(size(Tmain) ~= size(Tsub)) || any(size(Tmain) ~= size(Tsup))
    error 'Tmain, Tsup, Tsub are not identical in size'
  end
  if isempty(n) || length(n) > 1 ||  n < 1 || n ~= floor(n)
    error('n must be a positive scalar integer');
  end
  v = repmat(Tmain(:), n, 1);
  v = [v; repmat(Tsub(:), n - 1, 1)];
  v = [v; repmat(Tsup(:), n - 1, 1)];
  % now generate the index arrays. first the main diagonal
  [ind1, ind2, ind3] = ndgrid(0 : p - 1, 0 : q - 1, 0 : n - 1);
  rind = 1 + ind1(:) + p * ind3(:);
  cind = 1 + ind2(:) + q * ind3(:);
  % sub-diagonal
  [ind1, ind2, ind3] = ndgrid(0 : p - 1, 0 : q - 1, 0 : n - 2);
  rind = [rind; 1 + p + ind1(:) + p * ind3(:)];
  cind = [cind; 1 + ind2(:) + q * ind3(:)];
  % super-diagonal
  rind = [rind; 1 + ind1(:) + p * ind3(:)];
  cind = [cind; 1 + q + ind2(:) + q * ind3(:)];
  % build the final array all in one call to sparse
  T = sparse(rind, cind, v, n * p, n * q);
end
