function K_full = mrir_array_SMS_gre_recon(raw, K_IndivFull, kernel_size, ref_size, prot, PhaseShiftBtwSimulSlices, varargin)
%MRIR_ARRAY_SMS_GRE_RECON
%
% varargout = mrir_array_SMS_gre_recon(varargin)
%
%
% see also MRIR_ARRAY_SMS_GRE_KERNEL.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2013/dec/12
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  FLAG__fake_collapse = 0;
  if ( nargin >= 7 ),
    FLAG__fake_collapse = varargin{1};
  end;


  %==--------------------------------------------------------------------==%

  slicegroups = size(K_IndivFull{1}, 10);
  SMS = length(K_IndivFull);

  for slc = 1:slicegroups,

    disp(sprintf('slice %02d...', slc));
    K_Collapsed = 0;
    for band = 1:SMS,
      K_Indiv{band} = K_IndivFull{band}(:,:,:, 1,1, 1,1, 1,1,  slc);
    end;

    if ( FLAG__fake_collapse ),
      for count = 1:SMS,
        K_Collapsed = K_Collapsed+K_Indiv{count};
      end;
      K = MultisliceGRAPPA_1kernal_leakBlock(K_Collapsed,                    K_Indiv, kernel_size, ref_size, prot);
    else,
      K = MultisliceGRAPPA_1kernal_leakBlock(raw(:,:,:, 1,1, 1,1, 1,1, slc), K_Indiv, kernel_size, ref_size, prot);
    end;

    for SliceGroupCount = 1:SMS,
      K_full(:,:,:,:,:,:,:,:,:,slc+([SliceGroupCount-1]*slicegroups)) = CaipirinhaShift_K_v2(K(:,:,:,:,:,:,:,:,:,SliceGroupCount), SliceGroupCount,-PhaseShiftBtwSimulSlices);
    end;

  end;



  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
