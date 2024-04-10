function D = mrir_iDFT2_matrix(N1, N2, varargin)
%MRIR_IDFT2_MATRIX  efficient calculation of *unitary-scaled* 2D iDFT matrix
%
% D = mrir_iDFT2_matrix(N1, N2)
%
%
%  usage:
%
%     x = ifft2(k) * (sqrt(N1)*sqrt(N2));
%
%  is equivalent to
%
%     x = reshape(D*k(:), [N1, N2])

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/jan/11
% $Id: mrir_iDFT2_matrix.m,v 1.3 2009/10/26 03:27:52 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  FLAG__oversample_crop = 0;
  if ( nargin >= 3 ),
    FLAG__oversample_crop = varargin{1};
  end;


  %==--------------------------------------------------------------------==%

  % 'unitary' version of matrix, so adjoint(D) = hermitian(D)
  D1 = conj(dftmtx(N1))/sqrt(N1); D2 = conj(dftmtx(N2))/sqrt(N2);

  if ( exist('kronecker', 'file') ),
    disp('using Binary Singleton Expansion for calculation...');
    D = kronecker(D2, D1);
  else,
    D = kron(D2, D1);
  end;


  if ( FLAG__oversample_crop ),

    disp(sprintf('==> [%s]: removing rows from iDFT matrix to skip calculation of cropped voxels', mfilename));

    bnd = N1 / 4;

    ind_row = repmat( (bnd*1+1):(bnd*3-0), [N2,1]);
    ind_col = repmat([(0:(N2-1))*N1]', [1,N2]);
    ind = ind_row(:) + ind_col(:);

    D(ind,:) = [];

  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT2_matrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
