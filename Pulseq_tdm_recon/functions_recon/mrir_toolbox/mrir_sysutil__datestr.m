function [mydatestr, fmt] = mrir_sysutil__datestr(varargin)
%MRIR_SYSUTIL__DATESTR
%
% varargout = mrir_sysutil__datestr(varargin)
%
%
% see also DATESTR.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/jul/21
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  %if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  fmt = 'yyyy-mmm-dd HH:MM:SS PM';
  mydatestr = datestr(now, fmt);


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
