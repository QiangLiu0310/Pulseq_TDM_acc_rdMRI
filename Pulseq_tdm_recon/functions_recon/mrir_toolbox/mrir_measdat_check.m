function meas = mrir_measdat_check(meas)
%MRIR_MEASDAT_CHECK  check for idiosyncratic Siemens conventions and fix
%
% meas_fixed = mrir_measdat_check(meas)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/sep/20
% $Id: mrir_measdat_check.m,v 1.3 2012/05/25 09:53:22 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( ~isempty(regexp(meas.STATUS, 'checked')) ),
    warning('meas.dat file appears to already have been checked -- clear STATUS field to override');
    return;
  end;


  [meas, FLAG__fix_7T_channels] = mrir_hack__fix_7T_channels(meas);

  if ( FLAG__fix_7T_channels ),
    disp(sprintf('<i> [%s]: repaired meas.dat file---Bay 5 RCCS channel assignment issue with 32-channel coil FIXED', mfilename));
  end;

  msg = lastwarn;
  if ( strcmp(meas.STATUS, 'ABORTED') ),

    disp(sprintf('==> [%s]: repairing meas.dat file---discarding final repetition of aborted scan...', mfilename));

    %                                         1 2 3 4 5 6       7 8 9 0 1 2 3 4 5 6
    meas.data =           meas.data(          :,:,:,:,:,:,1:end-1,:,:,:,:,:,:,:,:,:);
    meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,1:end-1,:,:,:,:,:,:,:,:,:);

    if ( isfield(meas, 'evp') ),
      meas.evp.NRepMeas = mrir_ice_dimensions(meas.data, 'rep');
    end;

    disp(sprintf('<i> [%s]: truncation complete.', mfilename));
    meas.STATUS = 'TRUNCATED';

  end;


  if ( ~isfield(meas, 'prot') ),
    meas.STATUS = [meas.STATUS, '->checked'];
    return;
  end;


  NRep_read = meas.evp.NRepMeas;
  NRep_meas = NRep_read;
  meas.evp.NRepMeas = 1;

  % workaround for unusual Siemens convention #1:
  if ( meas.evp.NFirstRefLin == 0 & isfield(meas, 'patrefscan') ),
    meas.evp.NFirstRefLin = mrir_ice_dimensions(meas.patrefscan, 'lin') - meas.evp.NRefLin + 1;
    disp(sprintf('<i> [%s]: workaround for unusual Siemens convention applied (NFirstRefLin)', mfilename));
  end;

  if ( meas.evp.NFirstRefPar == 0 & isfield(meas, 'patrefscan') ),
    meas.evp.NFirstRefPar = mrir_ice_dimensions(meas.patrefscan, 'par') - meas.evp.NRefPar + 1;
    disp(sprintf('<i> [%s]: workaround for unusual Siemens convention applied (NFirstRefPar)', mfilename));
  end;


  % workaround for unusual Siemens convention #2:

  % (possibly Siemens fills in with GRAPPA fewer lines than are in FFT, so
  % last lines are effectively zero-padded; this could throw off SNR
  % calculations, so by overriding this we force "mrir_epi_GRAPPA" to fill
  % in same number of k-space lines as there are image lines.)
  if ( ~isempty(meas.evp.NAFLin) && (meas.evp.NAFLin == 1) && (prot.ucPhasePartialFourier == 1) && (meas.evp.NLinMeas < meas.evp.NImageLins) ),
    keyboard
    meas.evp.NLinMeas = meas.evp.NImageLins;
  end;

  meas.STATUS = [meas.STATUS, '->checked'];


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_measdat_check.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
