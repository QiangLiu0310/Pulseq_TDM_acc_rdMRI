function varargout = mrir_sysutil__depends(varargin)
%
% mrir_sysutil__depends(mfile)
%
% mrir_sysutil__depends(mfile1, mfile2, ...)
%
% toolbox_files = mrir_sysutil__depends(...)
%
% [toolbox_files, tar_commandline] = mrir_sysutil__depends(...)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/dec/28
% $Id: mrir_sysutil__depends.m,v 1.2 2009/12/28 20:34:24 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  p = path;
  addpath(fileparts(which(mfilename)));
  
  file_list = {};

  for ind = 1:nargin,

    mfile = varargin{ind};
%    deps = depfun(mfile, '-quiet');
    deps = matlab.codetools.requiredFilesAndProducts(mfile);

    file_list = unique({file_list{:}, deps{:}});

  end;

  [toolbox_file_mask, toolbox_priv_mask] = deal(zeros(1, length(file_list)));

  for ind = 1:length(file_list),

    pos = regexp(file_list{ind}, 'mrir_toolbox/');

    if ( ~isempty(pos) ),
      toolbox_file_mask(ind) = 1;
    end;

    pos = regexp(file_list{ind}, 'mris_toolbox/');

    if ( ~isempty(pos) ),
      toolbox_file_mask(ind) = 1;
    end;

    pos = regexp(file_list{ind}, 'jonp/matlab/');

    if ( ~isempty(pos) ),
      toolbox_file_mask(ind) = 1;
    end;

    pos = regexp(file_list{ind}, 'Contents');

    if ( ~isempty(pos) ),
      toolbox_file_mask(ind) = 0;
    end;

    pos = regexp(file_list{ind}, 'mrir_toolbox/private');

    if ( ~isempty(pos) ),
      toolbox_priv_mask(ind) = 1;
    end;

  end;



  toolbox_file_index = find(toolbox_file_mask);
  toolbox_priv_index = find(toolbox_priv_mask);

  % not used
  toolbox_call_index = setdiff(toolbox_file_index, toolbox_priv_index);

  % the files
  toolbox_files = strvcat(file_list{toolbox_file_index});

  
  %=---------------------------------------------------------------------==%
  % tar the file list
  
  tar_command_line = sprintf('%s ', file_list{toolbox_file_index});
  tar_command_line_relative = regexprep(tar_command_line, '\S*mrir_toolbox/', '');

  if ( length(nargin) >= 2 ),
    tarball = [tempname, datestr(now, '_yyyy_mm_dd')];
  else,
    tarball = [varargin{1}, datestr(now, '_yyyy_mm_dd')];
  end;

  tar_command_line_relative = sprintf('cd %s; tar Wcvphf %s/%s.tar %s; cd %s', fileparts(which(mfilename)), pwd, tarball, tar_command_line_relative, pwd);

  
  disp(toolbox_files);
  disp(tar_command_line_relative);

  
  mfile = varargin{1};
  
  mfile_tarball = sprintf('%s_%s', mfile, datestr(now, 30));
  
  fprintf(1, 'cd %s\n', pwd);
  fprintf(1, '\n\nmkdir -p %s\n\n', mfile);
  for ind = 1:length(toolbox_file_index),
    fprintf(1, 'cp -a %-120s %s/\n', file_list{toolbox_file_index(ind)}, mfile);
  end;
  fprintf(1, '\ntar Wcvlpf %s.tar %s/\n', mfile_tarball, mfile);
  
  
  
  if ( nargout >= 1 ),
    varargout{1} = toolbox_files;
    varargout{2} = tar_command_line_relative;
  end;

  path(p);

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__depends.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
