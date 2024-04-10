function [x, X] = mrir_iDFT2_covmtx(covmtx, imgsize, varargin)
%MRIR_IDFT2_COVMTX
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


  FLAG__method = 'fft';
  if ( nargin >= 3 ),
    FLAG__method = varargin{1};
  end;

  FLAG__crop = 1;
  if ( nargin >= 4 ),
    FLAG__crop = varargin{2};
  end;


  % number of k-space samples
  nK = unique(size(covmtx));

  % number of frequency-encode and phase-encode samples
  nF = imgsize(1); nP = imgsize(2);


  %==--------------------------------------------------------------------==%

  switch lower(FLAG__method),

   case {'econ', 1},
    %%% mrir_iDFT2_covariance__economy

    if ( FLAG__crop ),
      F = mrir_iDFT2_matrix(nF, nP, 1);
    else,
      F = mrir_iDFT2_matrix(nF, nP, 0);
    end;
      
    nI = size(F, 1);

    c = complex(zeros(nI,1), zeros(nI,1));

    for ind = 1:nI,

%      c(ind) = F(:,ind)' * covmtx * F(:,ind);
      c(ind) = conj(F(ind,:)) * covmtx * F(ind,:).';

    end;

    x = fftshift( reshape(c, [size(F,1)/nP, nP]) );

    % do not return full image-domain covariance when in economy mode
    X = [];

   case {'mtx_kron', 2},
    %%% mrir_iDFT2_covariance__dftmtx_kron

    F = mrir_iDFT2_matrix(nF, nP);

    X = F' * covmtx * F;
    x = fftshift( reshape(diag(X), [nF, nP]) );


   case {'mtx_loop', 3},
    %%% mrir_iDFT2_covariance__dftmtx_loop

    % "dftmtx.m" does not return unitary-scaled matrix, so need to scale here
    F = dftmtx2(nF, nP)/sqrt(nF)/sqrt(nP);
    X = F' * covmtx * F;
    x = fftshift( reshape(diag(X), [nF, nP]) );


   case {'fft', 'lilfft', 4},
    %%% mrir_iDFT2_covariance__fft

    for ind = 1:nK,

      %    if ( mod(ind, 100) == 0 ),
      %      disp(sprintf('%7d', ind));
      %    end;

      c = reshape(full(covmtx(:,ind)), imgsize);
      ic = mrir_iDFT(mrir_iDFT(c, 1), 2)/sqrt(nF)/sqrt(nP);

      R(:,ind) = sparse(ic(:));

    end;


    Rt = R';

    for ind = 1:nK,

      c = reshape(full(Rt(:,ind)), imgsize);
      ic = mrir_iDFT(mrir_iDFT(c, 1), 2)/sqrt(nF)/sqrt(nP);

      S(:,ind) = sparse(ic(:));

    end;

    X = S';

    x = reshape(diag(X), imgsize);

   case {'dotsum', 5},

    F = mrir_iDFT2_matrix(nF, nP, 1);
    x = sum( bsxfun(@times, F' * covmtx,  F.'), 2);

    X = [];

   case {'bigfft', 6},

    X = ...
        reshape(ifft(ifft(reshape(...
        reshape(ifft(ifft(reshape(...
            full(covmtx), ...
            [nF, nP, nF, nP]), [], 1), [], 2), [nK nK])', ...
            [nF, nP, nF, nP]), [], 1), [], 2), [nK nK])';

    x = fftshift( reshape(diag(X), [nF, nP]) );

   otherwise,

    error(sprintf('method "%s" not recognized', mat2str(FLAG__method)));

  end;  %% switch


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT2_covariance.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
