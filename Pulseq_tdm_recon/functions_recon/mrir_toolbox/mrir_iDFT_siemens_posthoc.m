function img_siemens = mrir_iDFT_siemens_posthoc(img_default)
%MRIR_IDFT_SIEMENS_POSTHOC  reformat data to match Siemens's iDFT results
%
% img_siemens = mrir_iDFT_siemens_posthoc(img_default)
%
% because siemens uses the fDFT (or "fft"), not the iDFT, for transforming
% k-space data into the image domain, their volumes come out flipped and
% shifted by one pixel in each Fourier Transformed dimension. this function
% flips and shifts the input data so that it matches the results of a
% Siemens reconstruction.
%
% this is intended for a post-hoc correction of image-domain data.
% alternatively, k-space data can be reconstructed the siemens way with
% "mrir_iDFT_siemens.m".
%
% see also MRIR_IDFT_SIEMENS.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/27
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % circshift and flip all Fourier Transformed dimensions
  % (if 2D dim(9) = 1 so this should have no effect)
  %         1  2 3 4 5 6 7 8  9 0 1 2 3 4 5 6
  shift = [-1,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0];

  img_siemens = flipdim(flipdim(flipdim(circshift(img_default, shift), 1), 2), 9);



  return;


  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
