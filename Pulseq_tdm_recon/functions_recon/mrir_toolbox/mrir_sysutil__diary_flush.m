function varargout = mrir_sysutil__diary_flush(varargin)
%MRIR_SYSUTIL__DIARY_FLUSH
%
% varargout = mrir_sysutil__diary_flush(varargin)
%
%
% see also MRIR_SYSUTIL__DIARY_FLUSH_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/dec/25
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  diary off; diary on;
  
  
  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
