function t0 = mrir_sysutil_alert
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/jun/27
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  t0 = now;
  
  beep; pause(0.5); beep; pause(0.5); beep; pause(0.1); beep

  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
