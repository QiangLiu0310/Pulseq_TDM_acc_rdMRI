function raw_grid = mrir_regrid(raw, prot)
%MRIR_REGRID
%
% varargout = mrir_regrid(varargin)
%
%
% see also MRIR_REGRID_TRAPEZOID.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/06
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  
  hyb_roft = mrir_iDFT_freqencode(raw);
  hyb_grid = mrir_regrid_trapezoid(hyb_roft, prot);
  raw_grid = mrir_fDFT_freqencode(hyb_grid);

  
  
  
  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
