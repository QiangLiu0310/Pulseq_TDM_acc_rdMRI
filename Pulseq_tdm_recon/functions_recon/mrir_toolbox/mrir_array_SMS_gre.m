function [gre, ref, raw] = mrir_array_SMS_gre(meas, varargin)
%MRIR_ARRAY_SMS_GRE
%
% varargout = mrir_array_SMS_gre(varargin)
%
%
% see also MRIR_ARRAY_SMS_EPI.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/dec/12
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  FLAG__phasestab = 1;
  if ( nargin >= 2 ),
    FLAG__phasestab = varargin{1};
    
    if ( ~isfield(meas, 'phasestabscan') ),
      warning('meas struct does not contain phasestabscan -- skipping phase stabilization');
      FLAG__phasestab = 0;
    end;
    
  end;
  
  FLAG__use_patrefscan = 1;
  if ( nargin >= 3 ),

    meas_sms1 = varargin{2};
    if ( ~isempty(meas_sms1) ),
      FLAG__use_patrefscan = 0;
    end;

  end;
  
  FLAG__deblur = 0;
  if ( nargin >= 4 ),
    FLAG__deblur = varargin{3};
  end;

  FLAG__ref_size_override = 0;
  if ( nargin >= 5 ),
    ref_size = varargin{4};
    
    FLAG__ref_size_override = 1;
  end;
  
  
  %==--------------------------------------------------------------------==%

  slcs = find(squeeze(meas.data(end/2,end/2,end, 1,1,1,1,1,1, :)));

  slicegroups = length(slcs);

  SMS = meas.prot.sWiPMemBlock_adFree(1);
  Nslices_after_recon = SMS * slicegroups;

  FOVshift = meas.prot.sWiPMemBlock_adFree(3);

  PhaseShiftBtwSimulSlices = 2*pi/FOVshift;

  SliceSep = meas.prot.dThickness/SMS;


  if ( FLAG__phasestab ),
    % mrir_array_SMS_gre_phasestab includes de-interleaving
    raw = mrir_array_SMS_gre_phasestab(meas, FOVshift);

  else,
    raw = mrir_image_slice_deinterleave(meas.data(:,:,:,1,1, 1,1, 1,1, slcs));
  end;


  NRefCol  = size(meas.patrefscan, 1);  %
  NRefLin  = size(meas.patrefscan, 2);  % = meas.evp.NRefLin;
  ref_size = [NRefCol, NRefLin];
  
  NImgCol  = size(meas.data, 1);  %
  NImgLin  = size(meas.data, 2);  % = meas.evp.NLinMeas
  
  if ( FLAG__use_patrefscan ),

    acs = mrir_zeropad(mrir_image_slice_deinterleave(meas.patrefscan), ...
                       [(NImgCol-NRefCol)/2,(NImgLin-NRefLin)/2,0, 0,0, 0,0, 0,0, 0], 'both');
  
  elseif ( FLAG__phasestab && isfield(meas_sms1, 'phasestabscan') ),
    acs = mrir_array_SMS_gre_phasestab(meas_sms1, FOVshift);
  else,
    acs = mrir_image_slice_deinterleave(meas_sms1.data);
  end;

  if ( prod(ref_size) == 1 ),
    ref_size = [size(acs, 1), size(acs, 2)];
  end;


  ref = mrir_array_SMS_ref(acs, SMS, slicegroups, PhaseShiftBtwSimulSlices);

  kernel_size = [3, 3, 2, 1];

  if ( FLAG__deblur ),
    disp(sprintf('==> [%s]: deblurring blipped data...', mfilename));
    raw = CaipirinhaDeblur_v2(raw, meas.prot, meas.evp, PhaseShiftBtwSimulSlices, SliceSep);
  end;
  
  gre = mrir_array_SMS_gre_recon(raw, ref, kernel_size, ref_size, meas.prot, PhaseShiftBtwSimulSlices);




  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
