% Load ScanArchive files and write to .h5 files for further processing

% get data file names and path
set_experimental_params;

% ghost calibration data
ifn = [datdir datfile_ghostcal]
ofn = [datdir ifn '.h5'];
if strcmp(scanner, 'GE')
    D = toppe.utils.loadsafile(ifn, 'acq_order', true);
else
    twix = mapVBVD(ifn);
    D = twix{2}.image.unsorted();
end
hmriutils.epi.io.draw2hdf(D, etl, np, ofn);

% slice GRAPPA calibration ('ACS') data
ifn = [datdir datfile_mb1]
ofn = [ifn '.h5'];
if strcmp(scanner, 'GE')
    % D = toppe.utils.loadpfile(ifn, [], [], [], 'acqOrder', true, 'returnAsDouble', false); 
    D = toppe.utils.loadsafile(ifn, 'acq_order', true);
else
    twix = mapVBVD(ifn);
    D = twix{2}.image.unsorted();
end
hmriutils.epi.io.draw2hdf(D, etl, np*mb, ofn);

% fMRI run(s)
ifn = [datdir datfile_rest]
ofn = [ifn '.h5'];
if strcmp(scanner, 'GE')
    D = toppe.utils.loadsafile(ifn, 'acq_order', true);
else
    twix = mapVBVD(ifn);
    D = twix{2}.image.unsorted();
end
% Only keep a few frames so the recon doesn't crash on a laptop
nfr = size(D,3)/etl/np;
nFramesToKeep = 4;
hmriutils.epi.io.draw2hdf(D(:,:,1:nFramesToKeep*etl*np), etl, np, ofn); %, 'maxFramesPerFile', 100);
