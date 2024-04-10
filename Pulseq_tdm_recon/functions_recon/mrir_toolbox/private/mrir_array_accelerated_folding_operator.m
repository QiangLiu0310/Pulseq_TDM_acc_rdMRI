function F = mrir_array_accelerated_aliasop(Nlin, R)
%MRIR_ARRAY_ACCELERATED_ALIASOP
%
% F = mrir_array_accelerated_aliasop(Nlin, R)

% jonathan polimeni <jonnyreb@padkeemao.nmr.mgh.harvard.edu>, 2007/may/04
% $Id: mrir_array_accelerated_folding_operator.m,v 1.1 2007/05/05 03:39:49 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % easy way to calculate aliasing operator: transform full-FOV image
  % indices to k-space, remove lines skipped during acceleration, then
  % transform back to image space.

  k_fullFOV = mrir_fDFT_freqencode(eye(Nlin));
  k_reduFOV = k_fullFOV(1:R:end, :);
  F = mrir_iDFT_freqencode(k_reduFOV)/Nlin*R;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/private/mrir_array_accelerated_folding_operator.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: