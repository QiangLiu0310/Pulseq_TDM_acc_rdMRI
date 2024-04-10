function [epi, sG, iG] = mrir_array_SMS_epi(meas, varargin)
%MRIR_ARRAY_SMS_EPI  reconstruct single repetiton of SMS-EPI data
%
% epi = mrir_array_SMS_epi(meas)
%
%   opt.ReadMultipleRepetitions = 0;
%   opt.SMSRefscan = 1;
%   meas = read_meas_dat(filename, opt);
%
% see also MRIR_ARRAY_SMS_KERNEL_RECON, MRIR_ARRAY_SMS_RECON_PARAMS, MRIR_EPI_GRAPPA_IMPROVED_DEV.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2015/jul/15
% $Id: mrir_array_SMS_epi.m,v 1.3 2015/10/17 00:32:17 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  FLAG__user_sG = 0;
  if ( nargin >= 2 && ~isempty(varargin{1}) ),
    FLAG__user_sG = 1
    sG = varargin{1};
  end;

  FLAG__user_iG = 0;
  if ( nargin >= 3 && ~isempty(varargin{2}) ),
    FLAG__user_iG = 1
    iG = varargin{2};
  end;

  FLAG__DO_LeakBlock = 1;
  if ( nargin >= 4 && ~isempty(varargin{3}) ),
    FLAG__DO_LeakBlock = varargin{3}
  end;


  %==--------------------------------------------------------------------==%
  %%% reconstruction options

  FLAG__LOCAL_PC = 1;

  % inplane-GRAPPA kernel size
  Nsrcx = 3;
  Nsrcy = 4;
  Nsrcz = 1;

  % inplane-GRAPPA regularization -- ON by default!
  eta = 1.0;
  chi = 0.01;

  sG_kernel_size = [5 5 2 1];
  FLAG__debug_with_fake_collapse = 0;


  %==--------------------------------------------------------------------==%

  meas = mrir_measdat_check(meas);

  NRep = mrir_ice_dimensions(meas.data, 'rep');


  [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(meas.prot, meas.evp);

  prot = meas.prot;
  prot.sSliceArray = meas.prot.sSliceArray(1:(NSlc/SMS));


  INTERLEAVE_REMOVE = 2;
  order = mrir_image_slice_deinterleave__sliceorder(NSlc, INTERLEAVE_REMOVE);

  if ( SMS == 1 ),
    epi = 0;
    return;
  end;


  %==--------------------------------------------------------------------==%

  if ( ~FLAG__user_sG ),
    smsrefscan         = mrir_array_SMS_recon_prune(meas.smsrefscan);
    if ( ~isequal(smsrefscan, meas.smsrefscan) && meas.prot.lAccelFactPE == 1 ),
      error('slice-GRAPPA ACS data preparation failed');
    end;

    smsrefscan_phascor = mrir_array_SMS_recon_prune(meas.smsrefscan_phascor);
    if ( ~isequal(smsrefscan_phascor, meas.smsrefscan_phascor) && meas.prot.lAccelFactPE == 1),
      error('slice-GRAPPA ACS navs preparation failed');
    end;

    [smsrefscan, check] = mrir_image_slice_deinterleave(smsrefscan);
    smsrefscan_phascor  = mrir_image_slice_deinterleave(smsrefscan_phascor);
    smsref_prep = mrir_evi_GRAPPA_prep__phascor_regrid(smsrefscan, smsrefscan_phascor, meas.prot);

    if ( ~isequal(order, check) ),
      error('slice order estimation failed');
    end;

    % add PhaseShift phase to the training data before collapsing (wrapper around "CaipirinhaShift_K_v2")
    smsref_phased = mrir_array_SMS_ref(smsref_prep, SMS, Ngroup, +PhaseShift);
  end;


  %==--------------------------------------------------------------------==%

  [dat_prep, patref_prep] = mrir_evi_GRAPPA_prep(meas, FLAG__LOCAL_PC); meas.data = []; meas.data_phascor1d = [];

  % remove phase shift from collapsed data (equivalent to setting phase of "center slice" to 0)
  dat_phased = CaipirinhaDeblur_v2_jrp(dat_prep, prot, meas.evp, PhaseShift, SliceSep); clear dat_prep


  %==--------------------------------------------------------------------==%

  for slc = 1:Ngroup,
    disp(sprintf('slice %02d...', slc));

    if ( FLAG__debug_with_fake_collapse ),
      dat_slc = 0;
      for count = 1:SMS,
        dat_slc = dat_slc + smsref_phased{count}(:,:,:, 1,1, 1,:, 1,1, slc);
      end;
    else,
      dat_slc = dat_phased(:,:,:, 1,1, 1,:, 1,1, slc);
    end;

    if ( ~FLAG__user_sG ),

      for band = 1:SMS,
        smsref{band} = smsref_phased{band}(:,:,:, 1,1, 1,1, 1,1,  slc);
      end;
      % wrapper around "MultisliceGRAPPA_*.m"
      [raw_phased, ~, sG{slc}] = mrir_array_SMS_kernel_recon(dat_slc,  smsref, sG_kernel_size, prot, FLAG__DO_LeakBlock);
    else,
      disp('using pre-computed slice-GRAPPA kernel!');
      raw_phased               = mrir_array_SMS_kernel_recon(dat_slc, sG{slc}, sG_kernel_size, prot, FLAG__DO_LeakBlock);
    end;

    % now subtract PhaseShift
    for band = 1:SMS,
      raw(:,:,:,:,:,:,:,:,:,slc+([band-1]*Ngroup)) = CaipirinhaShift_K_v2_jrp(raw_phased(:,:,:,:,:,:,:,:,:,band), band,-PhaseShift);
    end;
  end;


  %==--------------------------------------------------------------------==%

  if ( meas.prot.lAccelFactPE == 1 ),
    epi = mrir_conventional_2d(raw);
    iG = {};
    return;
  end;


  %==--------------------------------------------------------------------==%

  for slc = 1:Ngroup,

    acs = double(patref_prep(:,:,:,:,:,:,:,:,:,slc:Ngroup:end));

    for band = 1:SMS,

      % this is an abuse of notation -- acs_slc is either a volume or a kernelstruct
      if ( ~FLAG__user_iG ),
        acs_slc = acs(:,:,:,:,:,:,:,:,:,band);
      else,
        acs_slc = iG{slc};
      end;

      %                                               1 2 3 4 5 6  7 8 9  10
      [epi_rep, G] = mrir_epi_GRAPPA_improved_dev(raw(:,:,:,:,:,:,01,:,:,slc+([band-1]*Ngroup)), acs_slc, meas.prot, meas.evp, Nsrcx, Nsrcy, 1,1,1, eta, chi);

      %   1 2 3 4 5 6  7 8 9
      epi(:,:,:,:,:,:,01,:,:,slc+Ngroup*(band-1)) = epi_rep;

      iG{slc+Ngroup*(band-1)} = G;

      acs_img(:,:,:,:,:,:,1,:,:,slc+Ngroup*(band-1)) = mrir_epi_GRAPPA_improved_dev(smsref_prep(:,:,:,:,:,:,1,:,:,slc+([band-1]*Ngroup)), iG{slc+Ngroup*(band-1)}, meas.prot, meas.evp);

    end;  %% band
  end;  %% slc

  if ( NRep > 1 ),
    for rep = 2:NRep,
      disp(sprintf('==> [%s]: repetition %03d of %03d', mfilename, rep, NRep));
      for slc = 1:Ngroup,
        for band = 1:SMS,
          epi_rep = mrir_epi_GRAPPA_improved_dev(raw(:,:,:,:,:,:,rep,:,:,slc+([band-1]*Ngroup)), iG{slc+Ngroup*(band-1)}, meas.prot, meas.evp);
          epi(:,:,:,:,:,:,rep,:,:,slc+Ngroup*(band-1)) = epi_rep;
        end;  %% band
      end;  %% slc
    end;  %% rep
  end;



%  for slc = 1:NSlc,
%    acs_slc = double(patref_prep(:,:,:,:,:,:,:,:,:,slc));
%    [epi_rep, G{slc}] = mrir_epi_GRAPPA_improved_dev(smsref_cor(:,:,:,:,:,:,:,:,:,slc), ...
%                                                     acs_slc, meas.prot, meas.evp, ...
%                                                     Nsrcx, Nsrcy, 1,1,1, eta, chi);
%    epi(:,:,:,:,:,:,01,:,:,slc) = epi_rep;
%  end

  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SMS_epi.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
