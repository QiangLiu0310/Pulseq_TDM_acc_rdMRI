function img = mrir_array_SMS_epi(meas)
%MRIR_ARRAY_SMS_EPI
%
% varargout = mrir_array_SMS_epi(varargin)
%
%
% see also MRIR_ARRAY_SMS_EPI_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/01
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  % prep all repetitions
  [dat, acs] = mrir_evi_GRAPPA_prep(meas);
  
  raw = mrir_array_SMS_recon_LeakBlock(dat(:,:,:,1,1,1, 01, :, 1, :), acs, meas.prot, meas.evp);
  
  img = mrir_conventional_2d(raw);
  
  
  
  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
