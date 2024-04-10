function w = tukeywin(n,r)
%MRIR_MSTB_TUKEYWIN  copy of MATLAB's tukeywin.m (to avoid license useage)
%
% varargout = mrir_mstb_tukeywin(varargin)
%
%
% see also TUKEYWIN.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/nov/11
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  error(nargchk(1,2,nargin,'struct'));

  % default value for R parameter
  if ( nargin < 2 || isempty(r) ),
    r = 0.500;
  end;

  [n, w, trivialwin] = check_order(n);
  if ( trivialwin ), return, end;

  if ( r <= 0 ),
    w = ones(n,1);
  elseif ( r >= 1 ),
    w = hann(n);
  else,
    t = linspace(0,1,n)';
    % defines period of the taper as 1/2 period of a sine wave
    per = r/2;
    tl = floor(per*(n-1))+1;
    th = n-tl+1;
    % window is defined in three sections: taper, constant, taper
    w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
  end;


  return;

  
  
%**************************************************************************%
function [n_out, w, trivalwin] = check_order(n_in)

  w = [];
  trivalwin = 0;

  if ( ~(isnumeric(n_in) & isfinite(n_in)) ),
    error(generatemsgid('InvalidOrder'),'The order N must be finite.');
  end;

  % Special case of negative orders:
  if ( n_in < 0 ),
    error(generatemsgid('InvalidOrder'),'Order cannot be less than zero.');
  end;

  % Check if order is already an integer or empty
  % If not, round to nearest integer.
  if ( isempty(n_in) | n_in == floor(n_in) ),
    n_out = n_in;
  else,
    n_out = round(n_in);
    warning(generatemsgid('InvalidOrder'),'Rounding order to nearest integer.');
  end;

  % Special cases:
  if ( isempty(n_out) | n_out == 0 ),
    w = zeros(0,1);               % Empty matrix: 0-by-1
    trivalwin = 1;
  elseif ( n_out == 1 ),
    w = 1;
    trivalwin = 1;
  end;


  return;

  

  %************************************************************************%
  %%% $Source$
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
