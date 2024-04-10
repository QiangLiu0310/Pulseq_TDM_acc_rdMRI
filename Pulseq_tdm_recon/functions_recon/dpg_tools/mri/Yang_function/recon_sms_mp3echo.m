function [img_short, img_medium, img_long] = recon_sms_mp3echo( varargin )
%% add function:SEloop_mode 20200809
%  --------------------------------------------------------------------------------------------------
% |         SEloop_mode:0           |          SEloop_mode:1        |        SEloop_mode:1           | 
% |                                 |                               |                                |      
% | LongEcho | ShortEcho | Medium   | Medium | LongEcho | ShortEcho | ShortEcho | Medium  | LongEcho |  
% |  Part-I  |  Part-II  | Part-III | Part-I | Part-II  | Part-III  |    Part-I | Part-II | Part-III |
%  --------------------------------------------------------------------------------------------------
twix = varargin{1};
if iscell(twix)
    k = twix{length(twix)};
else
    k = twix;
end

%
SEloop_mode = varargin{2};   

%
tmp2  = varargin{3};
if iscell(tmp2)
    k_acs = tmp2{length(tmp2)};
else
    k_acs = tmp2;
end

%
NMultiplex = 3;% three echo sequence mp3;
hdr = parse_measdat_hdr(k.image.filename);
hdr_acs = parse_measdat_hdr(k_acs.image.filename);
if ((hdr.R~=hdr_acs.R)||(hdr.NRefLin~=hdr_acs.NRefLin)||(NMultiplex*hdr.NSlc~=hdr_acs.NSlc)||(hdr.SMS~=hdr_acs.SMS))
    error( '---------The dimension of externalACS data does not match------------' );
end

%
Ncoils = hdr.ncoils;
NRefLin = hdr.NRefLin;
NSlc = hdr.NSlc;
NRep = k.image.NRep;
NImgLin= length(k.image.Lin(1):hdr.R:k.image.NLin);
NLocPhz=3;%  3 lines at ky=0 are acquired, RO+ RO- RO+;
Ngroup=hdr.NSlc/hdr.SMS;
PhaseShift=2*pi/hdr.FOVshift;

slic_index_acs = hdr_acs.slcindx;
[~,order_acs]=sort(slic_index_acs);


%%%%%%%Caution here%%%%%%%%%%%%%%%%%%% determine the imaging order of MB group
starP=NImgLin*NMultiplex*NSlc+1;
endP =NImgLin*NMultiplex*NSlc+1+NImgLin*NMultiplex*(Ngroup-1);
intervalP=NImgLin*NMultiplex;
slicePos =k.image.slicePos(1:3,starP:intervalP:endP);
[~,order_MB]=sort(sum(slicePos,1));
%%%%%%%Caution here%%%%%%%%%%%%%%%%%%% determine the imaging order of MB group


% prot = eval_ascconv( k.image.filename );%Qiang
load('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/data_processing/data_processing_01_07/gslider_prot.mat')
evp.NSlcMeas = length(hdr.slcindx);

%% % for higher efficiency, high-spatital frequency EPI data is sampled 'on the
% ramps', meaning before the gradient trajectory is stable.  Thus, it
% needs to be 'gridded' onto a constant-velocity-gradient grid prior to
% processing. The VRGF approach does this, e.g. converting 256 sampled-grid
% points into 128 Cartesian-grid points.
%
v = extract_vrgf(k.image.filename);         % build a VRGF conversion
v_acs = extract_vrgf(k_acs.image.filename); 

NFreq_inres=k.phasecor.sqzSize(1);
NFreq_outres=size(v,2);


if hdr.R>2
    NSeg=hdr.R; %reference data to be acquired segmented
else
    NSeg=1;
end

sz_refscanPC=[NFreq_outres  Ncoils  NLocPhz          NMultiplex*NSlc   NSeg];
sz_refscan  =[NFreq_outres  Ncoils  NRefLin/NSeg     NMultiplex*NSlc   NSeg];

sz_phasecor =[NFreq_outres  Ncoils  NLocPhz        NMultiplex   Ngroup   NRep];
sz_image =   [NFreq_outres  Ncoils  NImgLin        NMultiplex   Ngroup   NRep];



%% GRAPPA calibration data extraction 

%refscanPC: ghost-correction data for k.refscan
refscanPC_raw_tmp  = k_acs.refscanPC.unsorted();
refscanPC_raw = reshape(v_acs'*reshape(refscanPC_raw_tmp,[NFreq_inres Ncoils*NLocPhz*NMultiplex*NSlc*NSeg]), sz_refscanPC);
% refscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
refscanPC_ROp_ROn(:,:,1,:,:,:) = squeeze((refscanPC_raw(:,:,1,:,:,:)+refscanPC_raw(:,:,3,:,:,:))/2);%ROp
refscanPC_ROp_ROn(:,:,2,:,:,:) = squeeze(refscanPC_raw(:,:,2,:,:,:));%ROn


%refscan: ACS data that can be used to train GRAPPA coefficients
refscan_raw_tmp = k_acs.refscan.unsorted();
refscan_raw = reshape(v_acs'*reshape(refscan_raw_tmp,[NFreq_inres Ncoils*NRefLin*NMultiplex*NSlc]), sz_refscan);

refscanPC  = squeeze(refscanPC_ROp_ROn(:,:,:,order_acs,:));
refscan    = squeeze(refscan_raw(:,:,:,order_acs,:));

refscan_pc_tmp     = PhaseCorrect_yang(refscanPC,refscan);
for i=1:NSeg
  refscan_pc(:,:,i:NSeg:NRefLin,:)=refscan_pc_tmp(:,:,:,:,i);  
end
C2=permute(refscan_pc,[1 3 4 2]);

C2_P1=C2(:,:,1+0*NSlc:1*NSlc,:);
C2_P2=C2(:,:,1+1*NSlc:2*NSlc,:);
C2_P3=C2(:,:,1+2*NSlc:3*NSlc,:);


k_trgt_P1 = C2_P1;
k_trgt_P2 = C2_P2;
k_trgt_P3 = C2_P3;

sz = size(k_trgt_P1);

if ~exist('slices','var'); slices = 1:hdr.NSlc; end;
if ~exist('chi','var'); chi = 1e-6; end;
if ~exist('eta','var'); eta = 1; end;


phzabs = exp(j*(PhaseShift*(0:(size(k_trgt_P1,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );

for cnt=0:(hdr.SMS-1),
  k_trgt_P1(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt_P1(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
  k_trgt_P2(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt_P2(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
  k_trgt_P3(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt_P3(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
end;



if ~exist('Npppp','var')
    
    fprintf('%d',mod(1:length(slices),10));fprintf('\n');
    for slc = 1:NSlc;
      fprintf('o');
      [~,~,Np_P1{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt_P1(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R*[1 2],...
                                       'chi',chi,'eta',eta );
                                   
      [~,~,Np_P2{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt_P2(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R*[1 2],...
                                       'chi',chi,'eta',eta );
                                   
      [~,~,Np_P3{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt_P3(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R*[1 2],...
                                       'chi',chi,'eta',eta );
    end;
    fprintf('\n');
 
end

z = [];
for cnt=1:length(Np_P1{slices(1)})
  dk = diff(find( Np_P1{slices(1)}(cnt).pattern == '*' ));
  z(dk) =  1;
end;

if ( (length(z)< 2*hdr.R) || ( z( 2*hdr.R ) == 0 ) )
  disp(['for this script to run properly, GRAPPA parameters in Np need to cover both R and 2R accelerations']);
  error(['need to regenerate Np parameters'])
end;

%% image data extraction
phasecor_raw_tmp   = k.phasecor.unsorted();
phasecor_raw =phasecor_raw_tmp(:,:,NLocPhz*NMultiplex*NSlc+1:end);
vrgf_phasecor =reshape(v'*reshape(phasecor_raw,[NFreq_inres Ncoils*NLocPhz*NMultiplex*Ngroup*NRep]), sz_phasecor);

% imgscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','3 echo','NSlc', 'NRep'
vrgf_phasecor_ROp_ROn(:,:,1,:,:,:)= squeeze((vrgf_phasecor(:,:,1,:,:,:)+vrgf_phasecor(:,:,3,:,:,:))/2);%ROp
vrgf_phasecor_ROp_ROn(:,:,2,:,:,:)= squeeze( vrgf_phasecor(:,:,2,:,:,:));%ROn
vrgf_phasecor_ROp_ROn=vrgf_phasecor_ROp_ROn(:,:,:,end:-1:1,:,:); %!!!!!local phase 321, image 123

%imgscan: ACS data that can be used to train GRAPPA coefficients
image_raw_tmp = k.image.unsorted();
image_raw=image_raw_tmp(:,:,NImgLin*NMultiplex*NSlc+1:end);
vrgf_image=reshape(v'*reshape(image_raw,[NFreq_inres Ncoils*NImgLin*NMultiplex*Ngroup*NRep]), sz_image);


Kimage= PhaseCorrect_yang(vrgf_phasecor_ROp_ROn,vrgf_image);

Kimage_short  =squeeze(Kimage(:,:,:,1,order_MB,:));% short
Kimage_medium =squeeze(Kimage(:,:,:,2,order_MB,:));% medium
Kimage_long   =squeeze(Kimage(:,:,:,3,order_MB,:));% long

%% recon
k_trgt_all{1}=k_trgt_P1;%region part1
k_trgt_all{2}=k_trgt_P2;%region part2
k_trgt_all{3}=k_trgt_P3;%region part3

Np_all{1}=  Np_P1;%region part1
Np_all{2}=  Np_P2;%region part2
Np_all{3}=  Np_P3;%region part3

for cmux=1:NMultiplex
    %--------------------------------------------------------------------
    index_mux = mod(SEloop_mode+cmux,NMultiplex)+1;
    k_trgt    = k_trgt_all{index_mux};
    Np        = Np_all{index_mux};

    %%%%calculate w, repeat Ngroup time
    for slc = 1:Ngroup
        iii = 1:size(k_trgt,2);
        sz = size(k_trgt);
        for cnt=1:hdr.SMS,
            in2{cnt} = zeros( sz([1 2 4]) );
            in2{cnt}(:,(sz(2)-length(iii))/2+iii,:) = squeeze( k_trgt(:,:,slc+(cnt-1)*Ngroup,:) );
        end;
        [~,w{slc}] =  MultisliceGRAPPA_2kernal_leakBlock( in2{1}, in2, [ 5 3 1 hdr.R ], 'full', prot);
        
    end
    %--------------------------------------------------------------------
    
    
    for cRep=1:NRep
               
        if cmux==1
            k_data_gc0 =permute(Kimage_short,[1 3 4 2 5]);
        elseif cmux==2
            k_data_gc0 =permute(Kimage_medium,[1 3 4 2 5]);
        else
            k_data_gc0 =permute(Kimage_long,[1 3 4 2 5]);
        end
        
        k_data_gc=k_data_gc0(:,:,:,:,cRep);
 
        
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
        for slc = 1:Ngroup;
            fprintf('.');
            
            in1 = zeros([ size(k_data_gc_deblur,1) nky*hdr.R sz(4)]);
            
            in1(:,hdr.R*(0:size(k_data_gc_deblur,2)-1)+1,:) = squeeze( k_data_gc_deblur(:,:,slc,:) );
            
            tmp = squeeze(MultisliceGRAPPA_2kernal_leakBlock( in1, w{slc}, [ 5 3 1 hdr.R ]) );
            
            
            %phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift) );
            %for cnt=1:hdr.SMS,
            %    tmp(:,:,:,cnt) = tmult( tmp(:,:,:,cnt), diag( phzabs.^(cnt-1)), 2);
            %end;
            
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
            
            %if verbose, keyboard; end;
        end
       fprintf('\n');
       k_all(:,:,:,:,cRep)=F2klb;
        
    end
     KKK{cmux}=k_all;
end


%% partial Fourier 
aa=KKK{1};
bb=KKK{2};
cc=KKK{3};
for cntsli=1:size(aa,3)
    
    for cntcoil=1:size(aa,4)
        for cRep=1:size(aa,5)
            img_aa(:,:,cntcoil,cRep)=flip(reconhd(aa(:,:,cntsli,cntcoil,cRep),size(aa,1),size(aa,2)));
            img_bb(:,:,cntcoil,cRep)=flip(reconhd(bb(:,:,cntsli,cntcoil,cRep),size(aa,1),size(aa,2)));
            img_cc(:,:,cntcoil,cRep)=flip(reconhd(cc(:,:,cntsli,cntcoil,cRep),size(aa,1),size(aa,2)));
            
            % flip should be done after reconhd, otherwise it will cause
            % very low SNR for large b-value, e.g. b=2000 or 3000
            %img_aa(:,:,cntcoil,cRep)=reconhd(flip(aa(:,:,cntsli,cntcoil,cRep)),size(aa,1),size(aa,2));
            %img_bb(:,:,cntcoil,cRep)=reconhd(flip(bb(:,:,cntsli,cntcoil,cRep)),size(bb,1),size(bb,2));
            %img_cc(:,:,cntcoil,cRep)=reconhd(flip(cc(:,:,cntsli,cntcoil,cRep)),size(cc,1),size(cc,2));
        end
        
    end
    img_short (:,:,cntsli,:) =sqrt(sum(abs( fif( img_aa ) ).^2,3));
    img_medium(:,:,cntsli,:) =sqrt(sum(abs( fif( img_bb ) ).^2,3));
    img_long  (:,:,cntsli,:) =sqrt(sum(abs( fif( img_cc ) ).^2,3));
end
