
% wrapper around Scott Hoge's NC-IGT Library implementation of Dual-Polarity
% GRAPPA EPI
%
% HOW TO USE:
%
%  A) enter into directory containing meas.dat files
%  B) make a subdirectory named after MID of meas.dat file to reconstruct
%     e.g., to recon "meas_MID101_ep2d_seg_se_1mm_iPAT20_seg5_TE55_*FID447.dat",
%           mkdir MID101; cd MID101
%  C) start matlab from MID directory


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>
% Tuesday, February 26, 2019  9:13:35 -0500
% Thursday, February 28, 2019 17:41:57 -0500
% Sunday, March 10, 2019 14:17:07 -0400

% based on: /cluster/visuo/users/scott/2019_02_21/MID236/runthis_jrp.m

% brief history:
% ver1) scott provides initial versions of "runthis.m", including:
%    /cluster/visuo/users/olivia/tmp/20190220_olivia_SEpilot014/MID101/runthis.m
%    /cluster/visuo/users/scott/2019_02_21/MID236/runthis.m
%
% ver2) jon adds support for proper NIfTI file creation, including (a)
% vox2ras information in header and (b) flipping/one-voxel-shifting to match
% Siemens DICOMs, and in the process discovers regridding bug and lack of
% apodization
%
% ver3) jon adds in regridding from the mrir_toolbox
%
% ver4) scott adds in support for optional GRAPPA wash, which also handles
% multi-FLEET data
%
% ver5) jon adds in user-controllable apodization from mrir_toolbox

% jonathan polimeni <jonp@arancine.nmr.mgh.harvard.edu>, 2019/mar/10
% $Id: mrir_example_script__EPI_recon_scott_DPG.m,v 1.2 2019/03/17 22:05:28 jonp Exp $
%**************************************************************************%

VERSION = '$Revision: 1.2 $';

addpath /cluster/visuo/users/share/dpg
dpgsetup


%==--------------------------------------------------------------------==%
%%% meas.dat file name

determine_fname;
if ~exist('k','var')
  k = mapVBVD(fname);
end;
v = extract_vrgf(fname);

[meas.prot, meas.evp] = read_meas_prot(fname);
regrid_trapezoid_prep = mrir_regrid_trapezoid_prep(meas.prot, size(k.refscan(:,:,:,1, 01,1,1,1,1,1,1), 1));

[~, measname] = fileparts(fname);


%==--------------------------------------------------------------------==%
%%% recon parameters

% set factor to 0.0 for no apodization
apodization_factor = 0.1;


hdr = parse_measdat_hdr(fname);
if ~isfield(hdr,'seg'), hdr.seg = 1; end;
hdr.R = hdr.R / hdr.seg ;
dpgkernel = sprintf('%dx5',2*hdr.seg);

% scott's "GRAPPA wash", in which kernel is used to reconstruct ACS data,
% then a second kernel is derived from this GRAPPA-reconstructed ACS data
hdr.grappa_wash = 1;

% normalization and regularization parameters
eta = 0.0;
chi = 1e-8;

savename = sprintf('%s__kern%s_eta%0.2g_chi%0.2g_apod0p%d_Gwash%d', ...
                   measname, dpgkernel, eta, chi, 10*apodization_factor, hdr.grappa_wash);


%==--------------------------------------------------------------------==%
%%% GRAPPA training

% k.refscan is ordered COL-CHA-LIN
% phzshift and recongrappa_multik expect LIN-COL-CHA

if ~exist('Ndpg','var')
for slc = 1:k.refscan.dataSize(5),
  fprintf('.');

  A = permute(mrir_regrid_trapezoid_scottdata(sum(k.refscan(:,:,:,1,slc,:,1,1,1,1,1),6), meas.prot, regrid_trapezoid_prep), [3,1,2]);

  B = permute(mrir_regrid_trapezoid_scottdata(sum(k.refscan(:,:,:,1,slc,:,1,1,1,1,2),6), meas.prot, regrid_trapezoid_prep), [3,1,2]);

  %A = permute(tmult( k.refscan(:,:,:,1,slc,1,1,1,1,1,1), v', 1),[3 1 2]);
  %B = permute(tmult( k.refscan(:,:,:,1,slc,1,1,1,1,1,2), v', 1),[3 1 2]);
  C = phzshift(A,B);

  if (hdr.grappa_wash),
    [~,~,Ng{slc}] = recongrappa_multik(size(C),C,vec(1:size(C,1)),'kernel','2x5','dks',2*hdr.R);

    for cnt=1:2*hdr.R,
      Fa(:,:,:,cnt) = recongrappa_multik(size(C),A(cnt:2*hdr.R:end,:,:),vec(cnt:2*hdr.R:size(C,1)), ...
                                         'kernel','2x5','dks',2*hdr.R,'N',Ng{slc});
      Fb(:,:,:,cnt) = recongrappa_multik(size(C),B(cnt:2*hdr.R:end,:,:),vec(cnt:2*hdr.R:size(C,1)), ...
                                         'kernel','2x5','dks',2*hdr.R,'N',Ng{slc});
      Fc(:,:,:,cnt) = recongrappa_multik(size(C),C(cnt:2*hdr.R:end,:,:),vec(cnt:2*hdr.R:size(C,1)), ...
                                         'kernel','2x5','dks',2*hdr.R,'N',Ng{slc});
    end;
    A2(:,:,slc,:) = sum(Fa,4)/size(Fa,4);
    B2(:,:,slc,:) = sum(Fb,4)/size(Fb,4);
    C2(:,:,slc,:) = sum(Fc,4)/size(Fc,4);
  else,
    A2(:,:,slc,:) = A;
    B2(:,:,slc,:) = B;
    C2(:,:,slc,:) = C;
  end;
  kin.p = squeeze(A2(:,:,slc,:));
  kin.n = squeeze(B2(:,:,slc,:));
  kin.target = squeeze(C2(:,:,slc,:));
  [~,~,Ndpg{slc}] = dpg_segepi(size(kin.p),kin,1:size(kin.p,1),'kernel',dpgkernel,'dks',hdr.R,'seg',hdr.seg,'normalize',eta,'chi',chi);

end;
fprintf('\n');

save_gparms(sprintf('%s__Ndpg.nd', savename), Ndpg);
end;

%==--------------------------------------------------------------------==%
%%% GRAPPA recon

for cntT = 1:hdr.seg:k.image.dataSize(9);

  for slc = 1:k.image.dataSize(5);
    % Fin.p and Fin.n are "zero padded" thus sparse and ordered
    %                                   r c p ?  s  ? ? ?   r                ? s

    Fin.p = permute(mrir_regrid_trapezoid_scottdata(sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,1),9), meas.prot, regrid_trapezoid_prep), [3,1,2]);

    Fin.n = permute(mrir_regrid_trapezoid_scottdata(sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,2),9), meas.prot, regrid_trapezoid_prep), [3,1,2]);


    %Fin.p = permute( tmult( sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,1),9), v', 1), [3 1 2]);
    %Fin.n = permute( tmult( sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,2),9), v', 1), [3 1 2]);

    Fdpg(:,:,slc,:) = dpg_recon( Fin, Ndpg{slc}, hdr.R, dpgkernel );

    if ( apodization_factor ),
      % apodize in RO and PE directions
      % Fdpga = mrir_filter_raw_apodize_1d(mrir_filter_raw_apodize_1d(Fdpg, 1, apodization_factor), 2, apodization_factor);
      % apodize in RO direction only
      Fdpga = mrir_filter_raw_apodize_1d(Fdpg, 1, apodization_factor);
    else,
      Fdpga = Fdpg;
    end;

  end

  % root-sum-of-squares combination
  Idpg(:,:,hdr.slcindx) = sqrt(sum(abs( fif(Fdpga) ).^2,4));
  %savend( sprintf('Idpg_%02d.nd',cntT), Idpg, 'flt' );

  vol = mrir_iDFT_siemens_posthoc(permute(Idpg, [2,1, 4,5,6,7,8,9,10, 3]));

  % save results as MGH file then convert to NIfTI using command-line tools
  % (e.g. mri_convert), since it is easier to generate an MGH header than a
  % NIfTI header in MATLAB
  mrir_save_mgh(sprintf('%s__vol%04d.mgh', savename, cntT), vol, meas.prot, meas.evp);

end;



%************************************************************************%
%%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_example_script__EPI_recon_scott_DPG.m,v $
%%% Local Variables:
%%% mode: Matlab
%%% fill-column: 76
%%% comment-column: 0
%%% End:
