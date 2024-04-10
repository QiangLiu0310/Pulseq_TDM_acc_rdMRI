function [dat_sms, order_sms] = mrir_array_SMS_recon_prune(dat, varargin)
%MRIR_ARRAY_SMS_RECON_PRUNE
%
% dat_sms = mrir_array_SMS_recon_prune(dat)
%
%
% see also MRIR_ARRAY_SMS_KERNEL_RECON

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/01
% $Id: mrir_array_SMS_recon_prune.m,v 1.1 2015/08/22 00:26:58 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( isempty(dat) ),
    dat_sms = dat;
    return;
  end;

  NSlc = mrir_ice_dimensions(dat, 'slc');
  
  dat_collapse = sum(sum(sum(sum(sum(sum(sum(sum(dat, 1), 3), 4), 5), 6), 7), 8), 16);

  dat_collapse(isnan(dat_collapse)) = 0;

  dat_lin = find(sum(dat_collapse, 10));
  dat_slc = find(sum(dat_collapse,  2));
  
  if ( isequal(vec(dat_slc), vec(1:NSlc)) ),
    warning('non-SMS data -- skipping slice pruning');
  end;

  dat_sms = dat(:,dat_lin,:,1,1,1,:, :, 1, dat_slc);


  if ( nargin >= 2 ),
    order = varargin{1};
    order_sms = order(dat_slc);
  else,
    order_sms = [];
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SMS_recon_prune.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
