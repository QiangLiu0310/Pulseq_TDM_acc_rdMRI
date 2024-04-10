function stop = mrir_sysutil__stopfile_check(varargin)
%MRIR_SYSUTIL__STOPFILE_CHECK
%
% stop = mrir_sysutil__stopfile_check()
% stop = mrir_sysutil__stopfile_check(stopfilename)
% stop = mrir_sysutil__stopfile_check(stopfilename, stopfilepath)
%
% default stopfile is "$HOME/MATLAB_stopfile"
%
%
% see also MRIR_SYSUTIL__TOUCHFILE.

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


  rehash
  stop = exist(stopfile, 'file');


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
