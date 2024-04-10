function result = mrir_sysutil__versioncheck(vercheck)
%MRIR_SYSUTIL__VERSIONCHECK
%
% results = mrir_sysutil__versioncheck(vercheck)
%
% a local copy of verLessThan, which only is shipped with versions more
% recent than 7.4.
%
% see also VERLESSTHAN.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/jul/26
% $Id: mrir_sysutil__versioncheck.m,v 1.1 2013/07/26 23:28:12 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  matlab_ver = sscanf(getfield(ver('matlab'), 'Version'), '%d.%d.%d')';
  matlab_ver(end+1:3) = 0;

  mrir_ver = sscanf(vercheck, '%d.%d.%d')';
  mrir_ver(end+1:3) = 0;

  result = (sign(matlab_ver - mrir_ver) * [1; .1; .01]) < 0;

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__versioncheck.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
