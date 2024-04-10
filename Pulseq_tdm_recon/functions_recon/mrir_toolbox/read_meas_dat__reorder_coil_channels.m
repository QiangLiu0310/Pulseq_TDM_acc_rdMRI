function varargout = read_meas_dat__reorder_coil_channels(header)
%READ_MEAS_DAT__REORDER_COIL_CHANNELS
%
%  [coil_index, coil_order] = read_meas_dat__reorder_coil_channels(header)

% to map using index:
%   dataval(coil_index);
% to map using order:
%   coil_order(dataind);


% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/27
% $Id: read_meas_dat__reorder_coil_channels.m,v 1.2 2012/05/26 09:58:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  % prepend "sCoilSelectMeas" to avoid confusion with "sCoilSelectUI" list.
  match = regexp(header, ['sCoilSelectMeas\S*\.sCoilElementID\.tCoilID\s*=\s*"(?<string>\w*)"'], 'names', 'once');
  coil_id = match.string;

  match = regexp(header, ['sCoilSelectMeas\.aRx\S*\.sCoilElementID\.tElement\s*=\s*"(?<string>\w*)"'], 'names');
  coil_element_matches = {match.string};
  [dummy, unique_index] = unique(coil_element_matches, 'first');
  coil_element_string = coil_element_matches(sort(unique_index));


  match = regexp(header, ['asCoilSelectMeas\S*\.lRxChannelConnected\s*=\s*(?<value>\w*)'], 'names');

  if ( isempty(match) ),
    match = regexp(header, ['sCoilSelectMeas\S*\.lADCChannelConnected\s*=\s*(?<value>\w*)'], 'names');
  end;

  coil_selected = str2num(str2mat({match(1:end/2).value}.'));


  % HACK: for some coils on 128-channel host, the banks of coil channels
  % are arranged such that the character in the middle of the string
  % should be sorted first, then the leading number, then the trailing
  % number. as a workaround, swap the leading number and the center
  % character so standard sort works appropriately.
  %  coil_element_mod = regexprep(coil_element_string, '([0-9]+)([a-zA-Z]+)(.*)', '$2$1$3');
  coil_element_mod = coil_element_string;

  [sorted, coil_index] = sort(coil_element_mod);
  [sorted, coil_order] = sort(coil_index);

  if ( nargout > 0 ),
    varargout{1} = coil_index;
    varargout{2} = coil_order;
    varargout{3} = coil_element_string;
  else,

    disp(sprintf('\t input channel  -->  sorted index'))
    for ind = 1:length(coil_index),
      disp(sprintf('\t     "%s"                %02d', ...
                   coil_element_string{ind}, coil_order(ind)));
    end;

  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_dat__reorder_coil_channels.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
