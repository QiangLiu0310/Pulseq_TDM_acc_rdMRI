function varargout = mrir_regrid_trapezoid_scottdata(raw, prot, varargin)
%MRIR_REGRID_TRAPEZOID_SCOTTDATA  wrapper for "mrir_regrid_trapezoid" to integrate into scott hoge's code
%
% raw_regrid = mrir_regrid_trapezoid_scottdata(raw, prot)
%
% parameters required: aflRegridADCDuration, alRegridRampupTime, alRegridRampdownTime

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2019/feb/28
% $Id: mrir_regrid_trapezoid_scottdata.m,v 1.1 2019/03/07 04:03:16 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  if ( nargin >= 3 ),
    prot_trapezoid = varargin{1};
  else,
    prot_trapezoid = [];
  end;

  DO__DEAPOD = 1;
  if ( nargin >= 4 ),
    DO__DEAPOD  = varargin{2};
  end;

  REGRID_NONE        = hex2dec('01');
  REGRID_TRAPEZOIDAL = hex2dec('02');
  REGRID_SINUSOIDAL  = hex2dec('04');


  %==--------------------------------------------------------------------==%

  % before becoming vulnerable to obscure error messages, make sure that
  % prot contains all the relevant parameter values
  if ( isempty(prot.aflRegridADCDuration) || ...
       isempty(prot.alRegridRampupTime)   || ...
       isempty(prot.alRegridRampdownTime) ),
    error('protocol is missing parameter values needed for regridding');
  end;


  % ASSUME: waveform is symmetric, i.e., rise time == fall time

  if ( isempty(prot_trapezoid) ),
    % calculate regriding parameters from sequence protocol parameter values
    prot_trapezoid = mrir_regrid_trapezoid_prep(prot, size(img, 1));
  end;

  if ( prot_trapezoid.alRegridMode == REGRID_NONE ),
    raw_regrid = raw;
    DO__DEAPOD = 0;
  else,
    raw_regrid = mrir_regrid_trapezoid_apply(raw, prot_trapezoid);
  end;

  clear raw;

  % project data back into x-space for deapodization AND cropping
  img_regrid = mrir_iDFT_freqencode(raw_regrid);


  if ( DO__DEAPOD ),
    img_deapod = mrir_regrid_trapezoid_rolloff(img_regrid, prot_trapezoid);

    img_final = mrir_image_crop(img_deapod, prot);

  else,
    %    disp('skipping deapodization');
    img_final = mrir_image_crop(img_regrid, prot);

  end;

  raw_final = mrir_fDFT_freqencode(img_final);

  varargout{1} = raw_final;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_regrid_trapezoid_scottdata.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
