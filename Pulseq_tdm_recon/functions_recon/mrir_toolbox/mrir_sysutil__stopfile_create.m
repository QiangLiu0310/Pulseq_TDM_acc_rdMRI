function stop = mrir_sysutil__stopfile_create(varargin)
%MRIR_SYSUTIL__STOPFILE_CREATE
%
% stop = mrir_sysutil__stopfile_create()
% stop = mrir_sysutil__stopfile_create(stopfilename)
% stop = mrir_sysutil__stopfile_create(stopfilename, stopfilepath)
%
% default stopfile is "$HOME/MATLAB_stopfile"
%
%
% see also MRIR_SYSUTIL__STOPFILE_CHECK

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2014/feb/22
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  %if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  stopfilename = 'MATLAB_stopfile';
  stopfile = fullfile(getenv('HOME'), stopfilename);

  if ( nargin >= 1 ),
    stopfile = varargin{1};
  end;

  if ( nargin >= 2 ),
    stopfile = fullfile(varargin{2}, varargin{1});
  end;

  fp = mrir_sysutil__sessionfile(stopfile);
  rehash


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
