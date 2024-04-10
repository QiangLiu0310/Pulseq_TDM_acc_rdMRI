function [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot, evp)
%MRIR_ARRAY_SMS_RECON_PARAMS
%
% [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot, evp)
%
%
% see also MRIR_ARRAY_SMS_EPI, MRIR_ARRAY_SMS_KERNEL_RECON.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/nov/01
% $Id: mrir_array_SMS_recon_params.m,v 1.2 2015/08/22 00:26:59 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if isfield(prot,'sSliceAcceleration')
    SMS = prot.sSliceAcceleration.lMultiBandFactor;
    FOVshift = prot.sSliceAcceleration.lFOVShiftFactor;
  else
  if isfield(prot,'sWiPMemBlock')
    SMS = prot.sWiPMemBlock.adFree(2); % number of slices in slice group
    FOVshift = prot.sWiPMemBlock.adFree(3); % the base FOV shift applied to the
                                          % SMS slices
  elseif isfield(prot,'sWipMemBlock')
    SMS      = prot.sWipMemBlock.adFree(5); % number of slices in slice group
    FOVshift = prot.sWipMemBlock.adFree(3); % the base FOV shift applied to the
                                            % SMS slices
  else
    error('WIP Mem Block is not defined');
  end;
  end;
  
  NSlc = evp.NSlcMeas; % number of measured slices
  Ngroup = NSlc/SMS;   % number of measured groups per volume

  if ( FOVshift == 1 ),
    PhaseShift = 0;
  else,
    PhaseShift = 2*pi/FOVshift;
  end;

  % 2015/oct/16
if nargout>5,
  for cnt=1:length( prot.sSliceArray.asSlice ),
    SlicePos(cnt,:) = [
        prot.sSliceArray.asSlice(cnt).sPosition.dSag;
        prot.sSliceArray.asSlice(cnt).sPosition.dCor;
        prot.sSliceArray.asSlice(cnt).sPosition.dTra;
                      ].';
  end;

  % the distance between the compressed slices:
  SliceSep = norm(SlicePos([NSlc/SMS]+1,:) - SlicePos(0+1,:));
end
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
