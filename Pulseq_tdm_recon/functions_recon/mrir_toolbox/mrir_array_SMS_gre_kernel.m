function K_IndivFull = mrir_array_SMS_gre_ref(K_IndivFull_nGrouped, SMS, slicegroups, PhaseShiftBtwSimulSlices)
%MRIR_ARRAY_SMS_GRE_REF
%
% K = mrir_array_SMS_gre_ref(acs, SMSfactor, Nslicegroups, PhaseShift)
%
%
% see also MRIR_ARRAY_SMS_GRE_RECON.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/dec/12
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Nslices_K_Collapsed = slicegroups;
  for SliceGroupCount = 1:SMS,
    K_IndivFull{SliceGroupCount} = CaipirinhaShift_K_v2(K_IndivFull_nGrouped(:,:,:,:,:,:,:,:,:,(1:Nslices_K_Collapsed)+(SliceGroupCount-1)*(Nslices_K_Collapsed) ),...
                                                  SliceGroupCount,PhaseShiftBtwSimulSlices);
  end;


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
