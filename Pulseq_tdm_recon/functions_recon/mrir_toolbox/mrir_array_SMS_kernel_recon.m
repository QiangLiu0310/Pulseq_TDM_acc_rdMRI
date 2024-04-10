function [raw, varargout] = mrir_array_SMS_kernel_recon(dat, ref, kernelsize, prot, varargin)
%MRIR_ARRAY_SMS_KERNEL_RECON  wrapper around "MultisliceGRAPPA"
%
% raw = mrir_array_SMS_kernel_recon(dat, ref, kernelsize, prot, varargin)
%
%
% see also MRIR_ARRAY_SMS_EPI.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2015/aug/13
% $Id: mrir_array_SMS_kernel_recon.m,v 1.1 2015/08/22 00:26:59 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  % todo: add options for 1-kernel vs. 2-kernel, regularization, etc.

  FLAG__DO_LeakBlock = 1;
  if ( nargin >= 5 && ~isempty(varargin{1}) ),
    FLAG__DO_LeakBlock = varargin{1};
  end;


  %==--------------------------------------------------------------------==%

  if ( FLAG__DO_LeakBlock ),
    if ( iscell(ref) ),
      % ACS data for training
      [raw, sG] = MultisliceGRAPPA_1kernal_leakBlock_jrp(dat, ref, kernelsize, 'full', prot);
    else,
      % pre-computed kernels
      [raw, sG] = MultisliceGRAPPA_1kernal_leakBlock_jrp(dat, ref, kernelsize);
    end;
  else,
    if ( iscell(ref) ),
      % ACS data for training
      [raw, sG] = MultisliceGRAPPA_1kernal_jrp(          dat, ref, kernelsize, 'full', prot);
    else,
      % pre-computed kernels
      [raw, sG] = MultisliceGRAPPA_1kernal_jrp(          dat, ref, kernelsize);
    end;
  end;


  if ( nargout >= 1 ),
    varargout{1} = raw;
    varargout{2} = sG;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SMS_kernel_recon.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
