function vol_flip = mrir_iDFT_siemens_posthoc_3d(vol)
%MRIR_IDFT_SIEMENS_POSTHOC_3D  wrapper around MRIR_IDFT_SIEMENS_POSTHOC for RxCxS volumes
%
% varargout = mrir_iDFT_siemens_posthoc_3d(varargin)
%
%
% see also MRIR_IDFT_SIEMENS_POSTHOC.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/apr/05
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  vol_flip = squeeze(mrir_iDFT_siemens_posthoc(permute(vol, [1,2,4,5,6,7,8,9,3])));

  
  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
