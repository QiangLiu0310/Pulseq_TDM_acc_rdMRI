function varargout = mrir_sysutil__touchfile(filestem, varargin)
%MRIR_SYSUTIL__TOUCHFILE
%
% touchfile = mrir_sysutil__touchfile(filestem, [ext], [msg], [FLAG__override])
%
%
% see also MRIR_SYSUTIL__WILDFILE.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/dec/12
% $Id: mrir_sysutil__touchfile.m,v 1.3 2013/06/16 19:45:15 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  ext = 'RUNNING';
  if ( nargin >= 2 ),
    ext = varargin{1};
  end;

  msg = '';
  if ( nargin >= 3 ),
    msg = varargin{2};
  end;

  FLAG__override = 0;
  if ( nargin >= 4 ),
    FLAG__override = varargin{3};
  end;


  global STOPIFERROR;  if ( isempty(STOPIFERROR) ),  STOPIFERROR = 1;   end;


  %==--------------------------------------------------------------------==%

  touchfile = sprintf('%s.%s', filestem, ext);

  if ( exist(touchfile, 'file') && ~FLAG__override ),
    wmsg = sprintf('WARNING: file "%s" exists; remove to proceed...', touchfile);
    disp(wmsg);
    enotify(wmsg);
    disp('    <<< dbcont to continue after removing file, return to abort >>>');
    mrir_sysutil__diary_flush(1);
    if ( STOPIFERROR ), varargout{1} = 0; keyboard; end;
  end;

  [status, result] = system(sprintf('touch %s', touchfile));
  if ( status ), error(result); end;

  fp = mrir_sysutil__sessionfile(touchfile, msg);


  varargout{1} = touchfile;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__touchfile.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
