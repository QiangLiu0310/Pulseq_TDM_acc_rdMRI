function varargout = mrir_sysutil__GB(variable)
%MRIR_SYSUTIL__GB  displays size of variable in GB -- without passing memory!
%
% GB = mrir_sysutil__GB(variable)
%
%
% see also WHOS, GB.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/oct/22
% $Id: mrir_sysutil__GB.m,v 1.1 2011/10/22 17:42:25 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  evalstr = sprintf('getfield(whos(''%s''), ''bytes'') / (2^30)', ...
                    inputname(1));


  gigabytes = evalin('caller', evalstr);

  space = repmat(' ', [1, 30 - length(inputname(1))]);

  if ( nargout == 0 ),
    disp(sprintf('  variable "%s":%s%10.4f gigabytes', inputname(1), space, gigabytes));
  else,
    varargout{1} = gigabytes;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__GB.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:



