function raw = mrir_array_SMS_recon_LeakBlock(dat, ref, prot, evp)
%MRIR_ARRAY_SMS_RECON_LEAKBLOCK
%
% raw = mrir_array_SMS_recon_LeakBlock(dat, ref, prot, evp)
%
%
% see also MRIR_ARRAY_SMS_RECON.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/oct/31
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  FLAG__fake_collapse = 0;


  addpath /autofs/space/bolt_001/users/kawin/MatlabCode/current_work/MultiSlice/Recon_Blip_v2 -end
  addpath(genpath('/autofs/space/bolt_001/users/kawin/MatlabCode/current_work/library/'), '-end')

  kernel_size = [5 5 2 1];


  %==--------------------------------------------------------------------==%

  [dat_sms, ref_sms, SMS, FOVshift, NSlc, Ngroup, PhaseShiftBtwSimulSlices] = mrir_array_SMS_recon_prep(dat, ref, prot, evp);


  %==--------------------------------------------------------------------==%


  for ind_slicegroup = 1:SMS,
    K_IndivFull{ind_slicegroup} = CaipirinhaShift_K_v2(ref_sms(:,:,:,:,:,:,:,:,:, ...
                                                  (1:Ngroup)+(ind_slicegroup-1)*Ngroup), ...
                                                  ind_slicegroup, ...
                                                  PhaseShiftBtwSimulSlices);
  end;


  %==--------------------------------------------------------------------==%


  for slc = 1:Ngroup,

    disp(sprintf('slice %02d...', slc));
    for band = 1:SMS,
      K_Indiv{band} = K_IndivFull{band}(:,:,:, 1,1, 1,1, 1,1,  slc);
    end;

    if ( FLAG__fake_collapse ),
      K_Collapsed = 0;
      for count = 1:SMS,
        K_Collapsed = K_Collapsed+K_Indiv{count};
      end;
      K = MultisliceGRAPPA_1kernal_leakBlock_jrp(K_Collapsed, K_Indiv, kernel_size, 'full', prot);
    else,
      K = MultisliceGRAPPA_1kernal_leakBlock_jrp(dat_sms(:,:,:, 1,1, 1,1, 1,1, slc), K_Indiv, kernel_size, 'full', prot);
    end;
    for SliceGroupCount = 1:SMS,
      raw(:,:,:,:,:,:,:,:,:,slc+([SliceGroupCount-1]*Ngroup)) = CaipirinhaShift_K_v2(K(:,:,:,:,:,:,:,:,:,SliceGroupCount), SliceGroupCount,-PhaseShiftBtwSimulSlices);
    end;

  end;

  img = mrir_conventional_2d(raw);


  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
