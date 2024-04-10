function filestr = mrir_sysutil__wildfile(wildfile, varargin)
%MRIR_SYSUTIL__WILDFILE
%
% filestr = mrir_sysutil__wildfile(wildstr, varargin)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/dec/28
% $Id: mrir_sysutil__wildfile.m,v 1.3 2012/01/19 16:15:29 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % expand wild cards
  [pathstr, basestr, ext] = fileparts(wildfile);

  % overwrite pathstr if input argument provided
  if ( nargin >= 2 ),
    pathstr = varargin{1};
  end;

  [status, listing] = unix( sprintf('ls -1 --color=never %s', fullfile(pathstr, [basestr, ext])) );
  
  if ( isempty(listing) ),
    disp(sprintf('pattern [%s] does not match existing meas file\n', wildfile));
    filestr = '';
    return;
  end;

  % take alphabetically first listing
  filestr = deblank(listing(1,:));

  
  
  
%  if ( ~exist(filestr, 'file') ),
%    disp(sprintf('pattern [%s] does not match existing file\n', wildfile));
%    filestr = '';
%    return;
%  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__wildfile.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


