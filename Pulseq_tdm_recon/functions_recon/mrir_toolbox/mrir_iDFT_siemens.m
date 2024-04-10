function ft = mrir_iDFT_siemens(raw, dim, varargin)
%MRIR_IDFT_SIEMENS  siemens-style inverse Discrete Fourier Transform
%
% ft = mrir_iDFT_siemens(raw, dim)
% ft = mrir_iDFT_siemens(raw, dim, N)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2012/nov/27
% $Id: mrir_iDFT.m,v 1.4 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  if ( nargin >= 3 ),
    Npoint = varargin{1};
  else,
    Npoint = size(raw, dim);
  end;


  %==--------------------------------------------------------------------==%

  % set FLAG__siemens_style to 1
  ft = mrir_iDFT(raw, dim, Npoint, 1);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
