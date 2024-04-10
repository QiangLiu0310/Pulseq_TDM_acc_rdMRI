function fp = mrir_sysutil__sessionfile(sessionfile, varargin)
%MRIR_SYSUTIL__SESSIONFILE
%
% varargout = mrir_sysutil__sessionfile(varargin)
%
%
% see also MRIR_SYSUTIL__SESSIONFILE_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2014/feb/22
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  msg = '';
  if ( nargin >= 2 ),
    msg = varargin{1};
  end;


  %==--------------------------------------------------------------------==%

  hostname = getenv('HOSTNAME');
  tty = getenv('TTY');

  [status, result] = system('/bin/ps | /bin/grep --color=never matlab_helper');
  if ( status ), error(result); end;
  ps = result;

  [fp, errormsg] = fopen(sessionfile, 'a');
  if ( fp == -1 ),
    error(errormsg);
  end;

  [diarypath, diaryfile, diaryext] = fileparts(get(0, 'DiaryFile'));
  if ( strcmp(get(0, 'Diary'), 'on') && ~isempty(diaryfile) ),
    msg = sprintf('%s\\n%s%s\\n', msg, diaryfile, diaryext);
  end;

  PBS_JOBID = getenv('PBS_JOBID');
  if ( ~isempty( PBS_JOBID ) ),
    msg = sprintf('%s\\n\\n[[ %s: %s, job "%s" ]]', msg, ...
                  getenv('PBS_SERVER'), getenv('PBS_JOBID'), getenv('PBS_JOBNAME'));
  end;

  fprintf(fp, '--------------------------------\n%s\n\n', mrir_sysutil__datestr());
  fprintf(fp, '%s\n%s\n%s\n', hostname, tty, ps);
  fprintf(fp, '%s\n\n\n', strrep(msg, '\n', char(10)));  % strrep replaces \n with newline character


  fclose(fp);


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
