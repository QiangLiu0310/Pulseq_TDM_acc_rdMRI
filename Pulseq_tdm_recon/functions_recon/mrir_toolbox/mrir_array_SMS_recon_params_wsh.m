function [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params_wsh(prot, evp)
%MRIR_ARRAY_SMS_RECON_PARAMS
%
% [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot, evp)
%
%
% see also MRIR_ARRAY_SMS_EPI, MRIR_ARRAY_SMS_KERNEL_RECON.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/01
% $Id: mrir_array_SMS_recon_params.m,v 1.2 2015/10/20 23:54:28 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  %==--------------------------------------------------------------------==%

  if ( ~isfield(prot, 'sWipMemBlock') || ...
       isempty(prot.sWipMemBlock.adFree) || ...
       length(prot.sWipMemBlock.adFree) < 3 ),
    
    SMS = 1;
    [FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = deal(0);
    return;
    
  end;

  
  SMS = prot.sWipMemBlock.adFree(3); if ( SMS == 0 ), SMS = 1; end;
  FOVshift = prot.sWipMemBlock.adFree(5);

  NSlc = evp.NSlcMeas;
  Ngroup = NSlc/SMS;

  if ( Ngroup < 1 || rem(NSlc, SMS) ),
    SMS = 1;
    [FOVshift, Ngroup, PhaseShift, SliceSep] = deal(0);
  end;
  
  
  if ( FOVshift == 1 ),
    PhaseShift = 0;
  else,
    PhaseShift = 2*pi/FOVshift;
  end;

  % 2015/oct/16
 for cnt=1:length( prot.sSliceArray.asSlice ),
    SlicePos(cnt,:) = [
        prot.sSliceArray.asSlice(cnt).sPosition.dSag;
        prot.sSliceArray.asSlice(cnt).sPosition.dCor;
        prot.sSliceArray.asSlice(cnt).sPosition.dTra;
                      ].';
  end;

  if ( SMS > 1 ),
    SliceSep = norm(SlicePos([NSlc/SMS]+1,:) - SlicePos(0+1,:));
  else,
    SliceSep = 0;
  end;

  % SliceSep = prot.dThickness/SMS;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SMS_recon_params.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

  
  
%        <ParamMap."sGroupArray"> 
%        {
%          
%          <ParamArray."asGroup"> 
%          {
%            <Default> <ParamMap.""> 
%            {
%              
%              <ParamDouble."distance"> 
%              {
%              }
%              
%              <ParamDouble."shift"> 
%              {
%              }
%              
%              <ParamDouble."thickness"> 
%              {
%              }
%            }
%            {  { 6.600  } { } { 59.400  } }   <--- skip factor (here 300%, 2 mm slice) & slice FOV
%            
%          }
%        }
