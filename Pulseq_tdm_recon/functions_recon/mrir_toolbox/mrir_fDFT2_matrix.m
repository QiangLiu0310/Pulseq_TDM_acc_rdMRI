function D = mrir_fDFT2_matrix(N1, N2)
%MRIR_FDFT2_MATRIX
%
% D = mrir_fDFT2_matrix(N1, N2)
%
%  
%  usage: 
%    
%     k = fft2(x);
%
%  is equivalent to
%
%     k = reshape(D*x(:), [N1, N2])

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/jan/11
% $Id: mrir_fDFT2_matrix.m,v 1.1 2009/01/12 00:47:19 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  
  D1 = dftmtx(N1); D2 = dftmtx(N2);
  D = kron(D2, D1);

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_fDFT2_matrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
