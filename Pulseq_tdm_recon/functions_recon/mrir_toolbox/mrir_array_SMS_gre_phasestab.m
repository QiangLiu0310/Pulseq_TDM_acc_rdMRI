function raw_stab = mrir_array_SMS_gre_phasestab(meas, FOVshift)
%MRIR_ARRAY_SMS_GRE_PHASESTAB
%
% varargout = mrir_array_SMS_gre_phasestab(varargin)
%
%
% see also MRIR_ARRAY_SMS_GRE.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/dec/12
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  slcs = find(squeeze(meas.data(end/2,end/2,end, 1,1,1,1,1,1, :)));

  raw_stab0 = zeros(size( meas.data ));

  if ( ~isfield(meas, 'phasestabscan') ),
    warning('data does not appear to have phase stabilization navigators -- skipping');
    raw_stab = mrir_image_slice_deinterleave(raw(:,:,:,1,1, 1,1, 1,1, slcs));
    return;
  end;
  
  disp(sprintf('==> [%s]: applying phase stabilization to FOV-shifted data', mfilename));
  
  for ind = 1:FOVshift,

    raw_stab0(:,ind:FOVshift:end,:, 1,1, 1,1, 1,1, :) = mrir_artifact_phase_stabilization(...
        meas.data(         :,ind:FOVshift:end,:, 1,1, 1,1, 1,1, :), ...
        meas.phasestabscan(:,ind:FOVshift:end,:, 1,1, 1,1, 1,1, :), ...
        meas.phasestabtime(:,ind:FOVshift:end,:, 1,1, 1,1, 1,1, :), ...
        meas.phasestabscan(:,ind,      :, 1,1, 1,1, 1,1, :), ...
        meas.prot.alTS*FOVshift, meas.prot.alTE);
  end;

  raw_stab = mrir_image_slice_deinterleave(raw_stab0(:,:,:,1,1, 1,1, 1,1, slcs));


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
