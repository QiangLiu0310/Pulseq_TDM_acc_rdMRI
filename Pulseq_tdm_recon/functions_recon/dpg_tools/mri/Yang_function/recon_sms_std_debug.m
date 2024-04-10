function [F2klb] = recon_sms_std_debug( varargin )
%% add function:SEloop_mode 20200809
twix_tmp  = varargin{1};
twix=twix_tmp{1,2};
hdr = parse_measdat_hdr(twix.image.filename);
Ncoils = hdr.ncoils;
NRefLin = hdr.NRefLin;
NSlc = hdr.NSlc;
R = hdr.R;
NEcho = 1;% std not mp;
NLocPhz=3;%  3 lines at ky=0 are acquired, RO+ RO- RO+;
NRep = twix.image.NRep;
slic_indexS = hdr.slcindx;
slic_indexM = hdr.slcindx;
slic_indexL = hdr.slcindx;

[~,reordr_S]=sort(slic_indexS);
[~,reordr_M]=sort(slic_indexM);
[~,reordr_L]=sort(slic_indexL);

slic_index_MB = [2:2:NSlc/hdr.SMS 1:2:NSlc/hdr.SMS];
[~,order_MB]=sort(slic_index_MB);

 prot = eval_ascconv( twix.image.filename );
evp.NSlcMeas = length(hdr.slcindx);

if isfield(hdr,'SMS')
  [SMS, FOVshift, NSlc, Ngroup, PhaseShift] = sms_recon_parms(hdr);
 
  [~,idx] = sort(hdr.slcindx);
  slc = 1:length(hdr.slcindx);
  indx = reshape( slc(idx), [ length(slc)/SMS SMS ])
  clear idx slc
 
  if (SMS>1)
    [~,idx] = sort(min(indx'));
  else
    [~,idx] = sort(indx);
  end;
  indx = indx(idx,:)
end;


if (nargin > 1)
   SEloop_mode = varargin{2};   
end

if (nargin > 2)
    twix_std_ = varargin{3};
    twix_std  = twix_std_{1,2};
    hdr_std = parse_measdat_hdr_Yang(twix_std.image.filename);
    if ((hdr.R~=hdr_std.R)||(hdr.NRefLin~=hdr_std.NRefLin)||(3*hdr.NSlc~=hdr_std.NSlc))
        disp( '---------The dimension of externalACS data does not match------------' );
        return;
    end
end


%% % for higher efficiency, high-spatital frequency EPI data is sampled 'on the
% ramps', meaning before the gradient trajectory is stable.  Thus, it
% needs to be 'gridded' onto a constant-velocity-gradient grid prior to
% processing. The VRGF approach does this, e.g. converting 256 sampled-grid
% points into 128 Cartesian-grid points.
%
v = extract_vrgf(twix.image.filename);         % build a VRGF conversion

NFreq_inres=twix.phasecor.sqzSize(1);
NFreq_outres=size(v,2);
NImgLin = size(twix.image.unsorted(),3)/(NEcho*(NSlc+Ngroup)*NRep);

if R>2
    NSeg=R; %reference data to be acquired segmented
else
    NSeg=1;
end

sz_refscanPC=[NFreq_outres  Ncoils  NLocPhz        NEcho   NSlc   NSeg];
sz_refscan  =[NFreq_outres  Ncoils  NRefLin/NSeg   NEcho   NSlc   NSeg];

sz_phasecor =[NFreq_outres  Ncoils  NLocPhz        NEcho   NSlc   NRep];
sz_image =   [NFreq_outres  Ncoils  NImgLin        NEcho   NSlc   NRep];



%% GRAPPA Reference data extraction 


%refscanPC: ghost-correction data for k.refscan
refscanPC_raw  = twix.refscanPC.unsorted();
vrgf_refscanPC = reshape(v'*reshape(refscanPC_raw,[NFreq_inres Ncoils*NLocPhz*NEcho*NSlc*NSeg]), sz_refscanPC);
% refscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
vrgf_refscanPC_ROp_ROn(:,:,1,:,:,:) = squeeze((vrgf_refscanPC(:,:,1,:,:,:)+vrgf_refscanPC(:,:,3,:,:,:))/2);%ROp
vrgf_refscanPC_ROp_ROn(:,:,2,:,:,:) = squeeze(vrgf_refscanPC(:,:,2,:,:,:));%ROn


%refscan: ACS data that can be used to train GRAPPA coefficients
vrgf_refscan_raw = twix.refscan.unsorted();
vrgf_refscan = reshape(v'*reshape(vrgf_refscan_raw,[NFreq_inres Ncoils*NRefLin*NEcho*NSlc]), sz_refscan);

vrgf_refscanPC_first   = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,reordr_S,:));
vrgf_refscan_first     = squeeze(vrgf_refscan(:,:,:,1,reordr_S,:));

refscan_first_temp     = PhaseCorrect_yang(vrgf_refscanPC_first,vrgf_refscan_first);
for i=1:NSeg
  refscan_first(:,:,i:NSeg:NRefLin,:)=refscan_first_temp(:,:,:,:,i);  
end
C2=permute(refscan_first,[1 3 4 2]);

if ~exist('verbose','var'); verbose = 0; end;
if ~exist('Ngroup','var'); Ngroup = hdr.NSlc / hdr.SMS; end;
if exist('C3')
  k_trgt = C3;
else
  k_trgt = C2;
end;
sz = size(k_trgt);

if ~exist('slices','var'); slices = 1:hdr.NSlc; end;
if ~exist('chi','var'); chi = 1e-6; end;
if ~exist('eta','var'); eta = 1; end;


phzabs = exp(j*(PhaseShift*(0:(size(k_trgt,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );

for cnt=0:(SMS-1),
  k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
end;


% collapsed data deinterleaving
if ( floor(hdr.NSlc/2) == (hdr.NSlc/2) )
  slcorder = [ 2:2:Ngroup 1:2:Ngroup ];
else
  slcorder = [ 1:2:Ngroup 2:2:Ngroup ];
end;

if ~exist('Np','var')
  if exist('Np.nd','file')
    Np = load_gparms('Np.nd');
  else
    fprintf('%d',mod(1:length(slices),10));fprintf('\n');
    for slc = slices; % 1:NSlc;
      fprintf('o');
      [~,~,Np{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R*[1 2],...
                                       'chi',chi,'eta',eta );
    end;
    fprintf('\n');
  end;
end;

z = [];
for cnt=1:length(Np{slices(1)})
  dk = diff(find( Np{slices(1)}(cnt).pattern == '*' ));
  z(dk) =  1;
end;

if ( (length(z)< 2*hdr.R) || ( z( 2*hdr.R ) == 0 ) )
  disp(['for this script to run properly, GRAPPA parameters in Np need to cover both R and 2R accelerations']);
  error(['need to regenerate Np parameters'])
end;

%% image data extraction
phasecor_raw   = twix.phasecor.unsorted();
sz_tmp=[NFreq_outres  Ncoils  NLocPhz     NEcho   (NSlc+Ngroup)   NRep];
phz_temp=reshape(v'*reshape(phasecor_raw,[NFreq_inres Ncoils*NLocPhz*NEcho*(NSlc+Ngroup)*NRep]), sz_tmp);
vrgf_phasecor  = phz_temp(:,:,:,:,NSlc+1:end);
% imgscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','3 echo','NSlc', 'NRep'
vrgf_phasecor_ROp_ROn(:,:,1,:,:)= squeeze((vrgf_phasecor(:,:,1,:,:)+vrgf_phasecor(:,:,3,:,:))/2);%ROp
vrgf_phasecor_ROp_ROn(:,:,2,:,:)= squeeze( vrgf_phasecor(:,:,2,:,:));%ROn

%imgscan: ACS data that can be used to train GRAPPA coefficients
vrgf_image_raw = twix.image.unsorted();
sz_tmp=[NFreq_outres  Ncoils  NImgLin     NEcho   (NSlc+Ngroup)   NRep];
img_tmp=reshape(v'*reshape(vrgf_image_raw,[NFreq_inres Ncoils*NImgLin*NEcho*(NSlc+Ngroup)*NRep]), sz_tmp);
vrgf_image = img_tmp(:,:,:,:,NSlc+1:end);


vrgf_phasecor_first = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,:,:));
vrgf_image_first    = squeeze(vrgf_image(:,:,:,:,:));

sz_k =[NFreq_outres  Ncoils  NImgLin     Ngroup   NRep];
Kimage_first_tmp = reshape(PhaseCorrect_yang(vrgf_phasecor_first,vrgf_image_first),sz_k);
Kimage_first =Kimage_first_tmp(:,:,:,order_MB);

if (twix.image.dataSize(9) > 1),
    Tcnt = 2;
else
    Tcnt = 1;
end;

%%-----------------------
% ksrc = twix;
% if ~exist('coils','var'); coils = 1:ksrc.image.dataSize(2); end;
% gg = tmult(ksrc.image(:,coils,:,1,:,:,1,1,Tcnt,:),v',1);
% vv = tmult( ksrc.phasecor(:,coils,:,1,:,:,1,1,Tcnt,:), v', 1 );
% 
% if (size(gg, length( size(gg) ))==2)
%     gg = sum( gg, length(size(gg)) );
% end;
% ggi = find( sum(sum( gg(:,:,:,1),1) ,2) );
% 
% sz = size(gg);
% sz = sz([1 3 5 2]);                   % [kx ky z coils]
% 
% v_data = permute( squeeze( vv(:,:,end,1,indx(:,1),:,:)),[1 4 3 2 5]);
% S1 = sum(v_data(:,:,:,:,1),2)./sum(v_data(:,:,:,:,1)~=0,2);
% S2 = sum(v_data(:,:,:,:,2),2)./sum(v_data(:,:,:,:,2)~=0,2);
% S1 = ffts(ifft(ffts( S1,1),[],1),1);
% S2 = ffts(ifft(ffts( S2,1),[],1),1);
% 
% 
% k_data = permute(gg(:,:,ggi,1,indx(:,1)),[1 3 5 2 4]);
% k_data(:,:,slcorder,:) = k_data;
% 
% 
% phz = comp_local_pc( S1, S2 );
% 
% k_data_gc_test = permute( phzapply( permute( k_data,[2 1 3 4]), phz ), [2 1 3 4] );
%%-----------------------
k_data_gc =permute(Kimage_first,[1 3 4 2]);
if isfield(prot,'sSliceAcceleration')
  % for data collected from the 'smsprod' seq, (eg 7T Terra data), data is already 'deblurred'
  k_data_gc_deblur = k_data_gc;
else
  % collapsed data Deblurring
  % this is needed for older 'mgh' versions of the sequence
  k_data_gc_deblur = CaipirinhaDeblur_v3_wsh( k_data_gc(:,:,:,:), prot, evp );
end;

if exist('runmod')
  if ( runmod == 1 )
    k_data_gc_deblur = k_data_gc_deblur(:,:,:,coils);
  end;
end;

nky = size(k_data_gc_deblur,2);

fprintf('%d',mod(1:length(slices),10)); fprintf('\n');
for slc = slices(1:end/hdr.SMS); % 1:Ngroup;
  fprintf('.');

  in1 = zeros([ size(k_data_gc_deblur,1) nky*hdr.R sz(4)]);

  in1(:,hdr.R*(0:size(k_data_gc_deblur,2)-1)+1,:) = squeeze( k_data_gc_deblur(:,:,slc,:) );

  if (~exist('w','var') || (length(w)<slc) || isempty(w{slc}) ),
    iii = 1:size(k_trgt,2);
    sz = size(k_trgt);
    for cnt=1:hdr.SMS,
      in2{cnt} = zeros( sz([1 2 4]) );
      in2{cnt}(:,(sz(2)-length(iii))/2+iii,:) = squeeze( k_trgt(:,:,slc+(cnt-1)*Ngroup,:) ); 
    end;
    [~,w{slc}] =  MultisliceGRAPPA_2kernal_leakBlock( in2{1}, in2, [ 5 3 1 hdr.R ], 'full', prot);
  end;

  tmp = squeeze(MultisliceGRAPPA_2kernal_leakBlock( in1, w{slc}, [ 5 3 1 hdr.R ]) );


 phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift) );
 for cnt=1:hdr.SMS,
   tmp(:,:,:,cnt) = tmult( tmp(:,:,:,cnt), diag( phzabs.^(cnt-1)), 2);
 end;

  phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );
  slcgrp = slc + [ 0:(hdr.SMS-1) ]*Ngroup;

  sz_in1 = size(in1);
  for cnt=1:hdr.SMS,
    curslc = slcgrp(cnt);
    Fa1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,1:2*hdr.R:end,:,cnt),[2 1 3]),1:2*hdr.R:sz_in1(2),'kernel','2x5','N',Np{curslc});
    Fb1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,(1+hdr.R):2*hdr.R:end,:,cnt),[2 1 3]),(1+hdr.R):2*hdr.R:sz_in1(2),'kernel','2x5','N',Np{curslc});
    
    Fa2 = permute(tunfold(fif(Fa1),2),[2 1]);
    Fb2 = permute(tunfold(fif(Fb1),2),[2 1]);

    Fb3 = phzshift( Fa2, Fb2,{'nofft','nocombo'} );
    Fb4 = ifi( trefold(permute(Fb3,[2 1]),sz_in1([2 1 3]),2) );

    Fc1 = zeros(size(Fa1));
    Fc1( 1:2*hdr.R:end, :, : ) = Fa1( 1:2*hdr.R:end, :, : );
    Fc1( (1+hdr.R):2*hdr.R:end, :, : ) = Fb4( (1+hdr.R):2*hdr.R:end, :, : );

    Fd1 = recongrappa_multik(size(Fc1),Fc1,[],'kernel','2x5','dks',hdr.R,'N',Np{curslc});

    F2klb(:,:,curslc,:) = Fd1;
    F2klb(:,:,curslc,:) =  tmult( F2klb(:,:,curslc,:), diag(conj(phzabs).^(cnt-1)), 1);
  end;

  if verbose, keyboard; end;
end;
fprintf('\n');
