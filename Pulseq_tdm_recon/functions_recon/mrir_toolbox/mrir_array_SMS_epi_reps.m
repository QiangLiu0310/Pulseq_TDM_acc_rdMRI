

opt.ReadMultipleRepetitions = 0;

meas = read_meas_dat(measfile, opt);



evp = meas.evp;
prot = meas.prot;

% extract position and acquisition parameters from protocol
mr_parms = mrir_measdat_parms(prot);
M0_vox2ras = mrir_measdat_vox2ras(prot, evp);

NRep = meas.evp.NRepMeas;


[epi_reps, sG, iG] = mrir_array_SMS_epi(meas);

opt.ReadMultipleRepetitions = 1;
opt.ExtractRepetition       = 1;

for rep = 2:NRep,

  disp(sprintf(' (repetition %03d of %03d)', rep, NRep));

  % read current repetition
  meas = read_meas_dat(measfile, opt, rep);
  meas.evp = evp;
  meas.prot = prot;

  epi = mrir_array_SMS_epi(meas, sG, iG);

  epi_reps = cat(7, epi_reps, epi);

end;



