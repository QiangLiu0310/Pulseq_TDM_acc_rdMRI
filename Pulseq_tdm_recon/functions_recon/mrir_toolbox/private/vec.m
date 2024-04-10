function m = vec(M)
%VEC  "vectorize" a matrix by converting it to the equivalent column vector
%
%
% m = vec(M)
%
% i.e., m == M(:) == reshape(M, '', 1)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/oct/10
% $Id: vec.m,v 1.1 2010/10/10 22:39:52 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  m = reshape(M, '', 1);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/jonp/MATLAB/vec.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


