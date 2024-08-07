function data_regrid = mrir_regrid_trapezoid_apply(data_orig, trapezoid)
%MRIR_REGRID_TRAPEZOID_APPLY  regrid ramp-sampled data with precomputed weights
%
% data_regrid = mrir_regrid_trapezoid_apply(data_orig, trapezoid)

% this code was poached directly from "ice_trapezoid_regrid.m", by Fa-Hsuan Lin.

% 2007/jan/03:

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/03
% $Id: mrir_regrid_trapezoid_apply.m,v 1.2 2011/03/18 01:59:45 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  REGRID_NONE        = hex2dec('01');
  REGRID_TRAPEZOIDAL = hex2dec('02');
  REGRID_SINUSOIDAL  = hex2dec('04');

  if ( trapezoid.alRegridMode == REGRID_NONE ),
    data_regrid = data_orig;
    return;
  end;

  if ( trapezoid.alRegridMode == REGRID_SINUSOIDAL ),
    error('sinusoidal regrid mode detected');
  end;

  
  %==--------------------------------------------------------------------==%

  %// Find the length of the vectors.
  lLenX = size(data_orig, 1);

  if ( lLenX ~= trapezoid.NxRegrid ),
    error('length of line for regridding does not match regrid function!');
  end;

  data_regrid = data_orig(:,:);

  for ind = 1:lLenX,

    jOrig = trapezoid.regrid_neighbor(ind,:)+1;

    data_regrid(ind,:) = trapezoid.regrid_convolve(ind,:)* ...
        (data_orig(jOrig,:)./ ...
         repmat(trapezoid.regrid_density(jOrig),[1,size(data_regrid,2)]));
  end;


  data_regrid = reshape(data_regrid, size(data_orig));


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_regrid_trapezoid_apply.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
