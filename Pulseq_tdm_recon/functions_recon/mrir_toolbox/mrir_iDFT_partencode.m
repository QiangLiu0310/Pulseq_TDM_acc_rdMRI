function paft = mrir_iDFT_partencode(raw, varargin)
%MRIR_IDFT_PARTENCODE  inverse Discrete Fourier Transform along partition encode
%
% paft = mrir_iDFT_partencode(raw);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id: mrir_iDFT_partencode.m,v 1.1 2007/09/17 20:18:36 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Npoint = size(raw, 9);
  if ( nargin >= 2 ),
    Npoint = varargin{2-1};
  end;

  if ( Npoint == 1 ),
    paft = raw;
    return;
  end;


  paft = mrir_iDFT(raw, 9, Npoint);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT_partencode.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
