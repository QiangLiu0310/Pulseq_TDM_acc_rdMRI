function raw_apod2 = mrir_filter_raw_apodize_2d(raw, varargin)
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/nov/05
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % choose a tukey tapering window [a.k.a. half-cosine taper] with
  % maximum unattenuated frequency at N/2*(1-k_taper_percent)
  k_taper_percent =  0.15;
  
  if ( nargin >= 3 ),
    k_taper_percent = varargin{1};
  end;
    
  
  raw_apod1 = mrir_filter_raw_apodize_1d(raw, mrir_DIM_COL, k_taper_percent);
  raw_apod2 = mrir_filter_raw_apodize_1d(raw, mrir_DIM_LIN, k_taper_percent);

  
  
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
