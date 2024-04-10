
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>
% Thursday, February 28, 2019 17:41:57 -0500

% /cluster/visuo/users/scott/2019_02_21/MID236/runthis_jrp.m


addpath /cluster/visuo/users/share/dpg
dpgsetup


determine_fname;
if ~exist('k','var')
  k = mapVBVD(fname);
end;
% v = extract_vrgf(fname);

[meas.prot, meas.evp] = read_meas_prot(fname);
regrid_trapezoid_prep = mrir_regrid_trapezoid_prep(meas.prot, size(k.refscan(:,:,:,1, 01,1,1,1,1,1,1), 1));



hdr = parse_measdat_hdr(fname);
if ~isfield(hdr,'seg'), hdr.seg = 1; end;
hdr.R = hdr.R / hdr.seg ;
dpgkernel = sprintf('%dx5',2*hdr.seg);

hdr.grappa_wash = 1;

% k.refscan is ordered COL-CHA-LIN
% phzshift and recongrappa_multik expect LIN-COL-CHA


for slc = 1:k.refscan.dataSize(5);
  fprintf('.');

  A = permute(mrir_regrid_trapezoid_scottdata(sum(k.refscan(:,:,:,1,slc,:,1,1,1,1,1),6), meas.prot, regrid_trapezoid_prep), [3,1,2]);

  B = permute(mrir_regrid_trapezoid_scottdata(sum(k.refscan(:,:,:,1,slc,:,1,1,1,1,2),6), meas.prot, regrid_trapezoid_prep), [3,1,2]);

  %A = permute(tmult( k.refscan(:,:,:,1,slc,1,1,1,1,1,1), v', 1),[3 1 2]);
  %B = permute(tmult( k.refscan(:,:,:,1,slc,1,1,1,1,1,2), v', 1),[3 1 2]);
  C = phzshift(A,B);
  
  if (hdr.grappa_wash)
    [~,~,Ng{slc}] = recongrappa_multik(size(C),C,vec(1:size(C,1)),'kernel','2x5','dks',2*hdr.R);
    
    for cnt=1:2*hdr.R;
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
  else
    A2(:,:,slc,:) = A;
    B2(:,:,slc,:) = B;
    C2(:,:,slc,:) = C;
  end;
  kin.p = squeeze(A2(:,:,slc,:));
  kin.n = squeeze(B2(:,:,slc,:));
  kin.target = squeeze(C2(:,:,slc,:));
  [~,~,Ndpg{slc}]=dpg_segepi(size(kin.p),kin,1:size(kin.p,1),'kernel',dpgkernel,'dks',hdr.R,'seg',hdr.seg,'normalize',0,'chi',1e-8);

end;
fprintf('\n');

save_gparms('Ndpg.nd',Ndpg);

for cntT = 1:hdr.seg:k.image.dataSize(9);

  for slc = 1:k.image.dataSize(5);
    % Fin.p and Fin.n are "zero padded" thus sparse and ordered
    %                                   r c p ?  s  ? ? ?   r                ? s

    Fin.p = permute(mrir_regrid_trapezoid_scottdata(sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,1),9), meas.prot, regrid_trapezoid_prep), [3,1,2]);

    Fin.n = permute(mrir_regrid_trapezoid_scottdata(sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,2),9), meas.prot, regrid_trapezoid_prep), [3,1,2]);


    %Fin.p = permute( tmult( sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,1),9), v', 1), [3 1 2]);
    %Fin.n = permute( tmult( sum(k.image(:,:,:,1,slc,1,1,1,cntT-1+(1:hdr.seg),1,2),9), v', 1), [3 1 2]);

    Fdpg(:,:,slc,:) = dpg_recon( Fin, Ndpg{slc}, hdr.R, dpgkernel );
  end

  Idpg(:,:,hdr.slcindx) = sqrt(sum(abs( fif(Fdpg) ).^2,4));
  savend( sprintf('Idpg_regrid_%02d.nd',cntT), Idpg, 'flt' );

  vol = mrir_iDFT_siemens_posthoc(permute(Idpg, [2,1, 4,5,6,7,8,9,10, 3]));

  mrir_save_mgh(sprintf('scott_recon_with_jon_regrid_%02d.mgh', cntT), vol, meas.prot, meas.evp);

end;
