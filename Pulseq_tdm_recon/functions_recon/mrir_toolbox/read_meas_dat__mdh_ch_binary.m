function mdh = read_meas_dat__mdh_ch_binary(binary)
%READ_MEAS_DAT__MDH_CH_BINARY
%  quickly parse VD channel header extracted as one binary vector
%
% mdh = read_meas_dat__mdh_ch_binary(binary)


% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense


% Anastasia Yendiki <ayendiki@nmr.mgh.harvard.edu>, 2011/oct/21
% $Id: read_meas_dat__mdh_ch_binary.m,v 1.2 2012/05/26 09:58:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  ind = 1;

  mdh.ulTypeAndChannelLength = typecastx(binary(ind:ind+3), 'uint32');
  ind = ind + 4;
  mdh.lMeasUID               = typecastx(binary(ind:ind+3), 'int32');
  ind = ind + 4;
  mdh.ulScanCounter          = typecastx(binary(ind:ind+3), 'uint32');
  ind = ind + 4;
  mdh.ulReserved1            = typecastx(binary(ind:ind+3), 'uint32');
  ind = ind + 4;
  mdh.ulSequenceTime         = typecastx(binary(ind:ind+3), 'uint32');
  ind = ind + 4;
  mdh.ulUnused2              = typecastx(binary(ind:ind+3), 'uint32');
  ind = ind + 4;
  mdh.ushChannelId           = typecastx(binary(ind:ind+1), 'uint16');
  ind = ind+2;
  mdh.ushUnused3             = typecastx(binary(ind:ind+1), 'uint16');
  ind = ind+2;
  mdh.ulCRC                  = typecastx(binary(ind:ind+3), 'uint32');
  ind = ind + 4;

  return;

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_dat__mdh_ch_binary.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

