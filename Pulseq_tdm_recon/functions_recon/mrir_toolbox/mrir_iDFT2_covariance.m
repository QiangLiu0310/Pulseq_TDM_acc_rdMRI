function [x, X] = mrir_iDFT2_covariance(covmtx, imgsize, FLAG__method)
%MRIR_IDFT2_COVARIANCE
%
% X = mrir_iDFT2_covariance(C, [m1, m2])
%
%   F = dftmtx(m1, m2)/sqrt(m1)/sqrt(m2);
%   c = F' * C * F
%   X = fftshift( reshape(diag(c), [m1, m2]) );

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/dec/11
% $Id: mrir_iDFT2_covariance.m,v 1.1 2008/12/12 02:56:07 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  disp('DEPRECATED! use "mrir_iDFT2_covmtx"');
  

  switch FLAG__method,

   case -1,
    %%% mrir_iDFT2_covariance__economy

    n = unique(size(covmtx));

    m1 = imgsize(1); m2 = imgsize(2);
    F = mrir_iDFT2_matrix(m1, m2);

    c = complex(zeros(n,1), zeros(n,1));
    
    for ind = 1:n,
      c(ind) = F(:,ind)' * covmtx * F(:,ind);
    end;

    x = fftshift( reshape(c, [m1, m2]) );
    X = [];

   case 0,
    %%% mrir_iDFT2_covariance__dftmtx_kron

    m1 = imgsize(1); m2 = imgsize(2);
    F = mrir_iDFT2_matrix(m1, m2);

    X = F' * covmtx * F;
    x = fftshift( reshape(diag(X), [m1, m2]) );

    
   case -2,
    %%% mrir_iDFT2_covariance__dftmtx_loop

    m1 = imgsize(1); m2 = imgsize(2);

    % "dftmtx.m" does not return unitary-scaled matrix, so need to scale here
    F = dftmtx2(m1, m2)/sqrt(m1)/sqrt(m2);
    X = F' * covmtx * F;
    x = fftshift( reshape(diag(X), [m1, m2]) );


   otherwise,
    %%% mrir_iDFT2_covariance__fft

    n = unique(size(covmtx));
    m1 = imgsize(1); m2 = imgsize(2);

    for ind = 1:n,

      %    if ( mod(ind, 100) == 0 ),
      %      disp(sprintf('%7d', ind));
      %    end;

      c = reshape(full(covmtx(:,ind)), imgsize);
      ic = mrir_iDFT(mrir_iDFT(c, 1), 2)/sqrt(m1)/sqrt(m2);

      R(:,ind) = sparse(ic(:));

    end;


    Rt = R';

    for ind = 1:n,

      c = reshape(full(Rt(:,ind)), imgsize);
      ic = mrir_iDFT(mrir_iDFT(c, 1), 2)/sqrt(m1)/sqrt(m2);

      S(:,ind) = sparse(ic(:));

    end;

    X = S';

    x = reshape(diag(X), imgsize);


  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT2_covariance.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
