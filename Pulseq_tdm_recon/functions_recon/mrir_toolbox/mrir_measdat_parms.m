function mr_parms = mrir_measdat_parms(prot)
%MRIR_MEASDAT_PARMS  extract parms from meas.prot for saving with save_mgh.m
%
% mr_parms = mrir_measdat_parms(prot)
%
%
% see also MRIR_MEASDAT_PROTOCOLNAME, MRIR_MEASDAT_VOX2RAS.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/27
% $Id: mrir_measdat_parms.m,v 1.1 2012/12/25 20:10:52 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  TR_ms = prot.alTR(1)/1000;
  TE_ms = prot.alTE(1)/1000;
  FA_rad = deg2rad(prot.adFlipAngleDegree(1));


  TI_ms = -1;
  if ( ~isempty(prot.alTI) ),
    TI_ms = prot.alTI(1)/1000;
    if ( prot.alTI == 300000 ),
      TI_ms = -1;
    end;
  end;
    
  mr_parms = [TR_ms, FA_rad, TE_ms, TI_ms];

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_measdat_parms.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
