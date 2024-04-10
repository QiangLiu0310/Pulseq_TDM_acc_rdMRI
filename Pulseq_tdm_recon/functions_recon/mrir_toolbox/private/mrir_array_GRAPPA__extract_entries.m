function [m, X, Y, C] = mrir_array_GRAPPA__extract_entries(M, x, y, c)


  [X, Y, C] = ndgrid(x,y,c);
  ind = sub2ind(size(M), X, Y, C);
  m = reshape(M(ind), [], 1);


  return;
