function [aliasmap, varargout] = mrir_array_accelerated_aliasmap(Nlin, R, ind_line)
%MRIR_ARRAY_ACCELERATED_ALIASMAP
%
% aliasmap = mrir_array_accelerated_aliasmap(Nlin, R, ind_line)

% jonathan polimeni <jonnyreb@padkeemao.nmr.mgh.harvard.edu>, 2007/may/04
% $Id: mrir_array_accelerated_aliasmap.m,v 1.1 2007/05/05 03:37:52 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % "aliasmap_full2redu" is a mapping from full-FOV pixel indices to
  % reduced-FOV pixel indices caused by subsampling during acceleration; and
  % "aliasmap" is the inverse mapping.

  % aliasmap_full2redu:      1 x Nlin  (each FULL maps to one REDU)
  % aliasmap          : Nlin/R x R     (each REDU maps to R FULL)

  % shift map by half a period if R is even
  if ( mod(R, 2) == 0 ),
    shift = floor(floor(Nlin/R) / 2);
  else,
    shift = 0;
  end;

  % the mod operator provides the aliasing mapping
  aliasmap_full2redu = mod( [1:Nlin] + shift, floor(Nlin/R));
  aliasmap_full2redu(find(aliasmap_full2redu==0)) = floor(Nlin/R);

  % invert table to map reduced FOV indices to full FOV indices

  % (since Nlin might not be evenly divisible by R, only find the first R matches)
  for ind = 1:(floor(Nlin/R)),
    aliasmap(ind, 1:R) = find(aliasmap_full2redu == ind, R);
  end;

  % calculate and return aliasop if requested
  if ( nargout >= 2 ),
    aliasop = zeros(R, Nlin);
    aliasop(sub2ind([R, Nlin], 1:R, aliasmap(ind_line,:))) = 1;
    varargout{1} = aliasop;
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/private/mrir_array_accelerated_aliasmap.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: