% Load Siemens .dat files and write to .h5 files for further processing

addpath ~/github/jfnielsen/scanLog/SiemensRead/

% get data file names and path
set_experimental_params;

% 3D GRE B0 map
if false
    ifn = [datdir datfile_b0];
    ofn = [ifn '.h5'];
    twix = mapVBVD(ifn);
    D = twix{2}.image.unsorted();          % [nx nc ny*(nz+1)]
    D = D(:,:,(2*params.b0.N(2)+1:end));   % discard dummy TRs (one kz loop)
    D = reshape(D, size(D,1), size(D,2), 2*params.b0.N(2), params.b0.N(3));  % [nx nc 2*ny nz]
    D = permute(D, [1 3 4 2]);             % [nx 2*ny nz nc]
    [~, mag] = toppe.utils.ift3(D(:,1:2:end,:,:));
end

% ghost calibration data (cal.dat)
ifn = [datdir datfile_cal];
ofn = [ifn '.h5'];
twix = mapVBVD(ifn);
D = twix{2}.image.unsorted();
hmriutils.epi.io.draw2hdf(D, etl, np, ofn);

% 2d reference scan (2d.dat) (used for slice GRAPPA)
ifn = [datdir datfile_mb1];
ofn = [ifn '.h5'];
twix = mapVBVD(ifn);
D = twix{2}.image.unsorted();
hmriutils.epi.io.draw2hdf(D, etl, np*mb, ofn);

% fMRI run(s)
ifn = [datdir datfile_mb6];
ofn = [ifn '.h5'];
twix = mapVBVD(ifn);
D = twix{2}.image.unsorted();
hmriutils.epi.io.draw2hdf(D, etl, np, ofn);
%nfr = 20;
%hmriutils.epi.io.draw2hdf(D(:,:,1:etl*np*nfr), etl, np, ofn);

