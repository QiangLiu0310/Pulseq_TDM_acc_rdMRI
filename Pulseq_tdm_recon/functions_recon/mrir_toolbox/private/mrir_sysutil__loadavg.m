function loadavg = mrir_sysutil__loadavg(varargin)
%MRIR_SYSUTIL__LOADAVG
%
% loadavg = mrir_sysutil__loadavg

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/15
% $Id: mrir_sysutil__loadavg.m,v 1.1 2008/12/16 01:51:24 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  %if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  filename = '/proc/loadavg';
  
  load(filename);
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/private/mrir_sysutil__loadavg.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
