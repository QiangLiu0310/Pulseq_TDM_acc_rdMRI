function filename = mrir_sysutil__screenshot(varargin)
%MRIR_SYSUTIL__SCREENSHOT  takes a full-screen screenshot using "import"
%
% pngfile = mrir_sysutil__screenshot(path)
% pngfile = mrir_sysutil__screenshot(filename)
%
%
% see also PRINTALL

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2014/jul/19
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin < 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % note: this is user:jonp specific
  SPACEpath = getenv('SPACE');
  MATLABDIRpath = [SPACEpath, '/MATLAB'];

  SCREENSHOTpath = [MATLABDIRpath, '/.SCREENSHOTS'];

  filename = '';
  if ( nargin >= 1 ),
    filename = varargin{1};

    if ( exist(filename, 'dir') ),
      userpath = filename;
      filename = '';
    else,
      userpath = fileparts(filename);
    end;

    if ( ~isempty(userpath) ),
      SCREENSHOTpath = userpath;
    end;
  end;

  if ( isempty(filename) ),

    logfilename = get(0, 'DiaryFile');
    TTY = regexprep(getenv('TTY'), '/', '-');

    filename = sprintf('screenshot__MATLAB_%s_%s_%s.png', ...
                       regexprep(getenv('HOSTNAME'), '\..*', ''), ...
                       datestr(now, 30), TTY);


  end;

  filename = fullfile(SCREENSHOTpath, filename);


  [status, result] = system(sprintf('import -window root %s', filename));
  if ( status ),
    error(result);
  end;


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
