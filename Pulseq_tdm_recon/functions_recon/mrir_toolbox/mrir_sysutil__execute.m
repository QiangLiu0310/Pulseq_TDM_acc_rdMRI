function varargout = mrir_sysutil__execute(instr)
%MRIR_SYSUTIL__EXECUTE
%
% varargout = mrir_sysutil__execute(varargin)
%
%
% see also SYSTEM.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/06
% $Id: mrir_sysutil__execute.m,v 1.1 2015/05/09 22:28:18 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [status, result] = system(instr);

  if ( status ),
    error(result);
  end;

  disp(result);
  

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__execute.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
