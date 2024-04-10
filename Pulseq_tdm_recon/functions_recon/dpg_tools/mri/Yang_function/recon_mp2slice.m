function [A] = recon_mp2slice( varargin )
%%
twix_  = varargin{1};
twix=twix_{1,2};
hdr = parse_measdat_hdr_Yang(twix.image.filename);
Ncoils = hdr.Ncoils;
NRefLin = hdr.NRefLin;
NSlc = hdr.NSlc;
R = hdr.R;
NEcho = 2;% two slice sequence mp2slice;
NPhaCorrLin=3;%  3 lines at ky=0 are acquired, RO+ RO- RO+;
NRep = twix.image.NRep;
% slic_indexS = [ 2 4 6 1 3 5];
% slic_indexL = [ 5 1 3 4 6 2];
slic_indexP1 = hdr.slcindx;
slic_indexP2 = hdr.slcindx;

[~,reordr_P1]=sort(slic_indexP1);
[~,reordr_P2]=sort(slic_indexP2);


if (nargin > 1)
    twix_acs_ = varargin{2};
    twix_acs  = twix_acs_{1,2};
    hdr_acs = parse_measdat_hdr_Yang(twix_acs.image.filename);
    if ((hdr.R~=hdr_acs.R)||(hdr.NRefLin~=hdr_acs.NRefLin)||(2*hdr.NSlc~=hdr_acs.NSlc))
        disp( 'The dimension of externalACS data does not match' );
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
v_acs = extract_vrgf(twix_acs.image.filename);

NFreq_inres=twix.phasecor.sqzSize(1);
NFreq_outres=size(v,2);
NImgLin = size(twix.image.unsorted(),3)/(NEcho*NSlc*NRep);

if R>2
    NSeg=R; %reference data to be acquired segmented
else
    NSeg=1;
end

sz_refscanPC=[NFreq_outres  Ncoils  NPhaCorrLin    NEcho   NSlc   NSeg];
sz_refscan  =[NFreq_outres  Ncoils  NRefLin/NSeg   NEcho   NSlc   NSeg];

sz_phasecor =[NFreq_outres  Ncoils  NPhaCorrLin    NEcho   NSlc   NRep];
sz_image =   [NFreq_outres  Ncoils  NImgLin        NEcho   NSlc   NRep];



%% GRAPPA Reference data extraction 

    
%refscanPC: ghost-correction data for k.refscan
refscanPC_raw    = twix.refscanPC.unsorted();
vrgf_refscanPC   = reshape(v'*reshape(refscanPC_raw,[NFreq_inres Ncoils*NPhaCorrLin*NEcho*NSlc*NSeg]), sz_refscanPC);
% refscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
vrgf_refscanPC_ROp_ROn(:,:,1,:,:,:)= squeeze((vrgf_refscanPC(:,:,1,:,:,:)+vrgf_refscanPC(:,:,3,:,:,:))/2);%ROp
vrgf_refscanPC_ROp_ROn(:,:,2,:,:,:)= squeeze(vrgf_refscanPC(:,:,2,:,:,:));%ROn


%refscan: ACS data that can be used to train GRAPPA coefficients
vrgf_refscan_raw = twix.refscan.unsorted();
vrgf_refscan     = reshape(v'*reshape(vrgf_refscan_raw,[NFreq_inres Ncoils*NRefLin*NEcho*NSlc]), sz_refscan);


%% image data extraction
phasecor_raw   = twix.phasecor.unsorted();
vrgf_phasecor  = reshape(v'*reshape(phasecor_raw,[NFreq_inres Ncoils*NPhaCorrLin*NEcho*NSlc*NRep]), sz_phasecor);
% imgscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
vrgf_phasecor_ROp_ROn(:,:,1,:,:,:) = squeeze((vrgf_phasecor(:,:,1,:,:,:)+vrgf_phasecor(:,:,3,:,:,:))/2);%ROp
vrgf_phasecor_ROp_ROn(:,:,2,:,:,:) = squeeze(vrgf_phasecor(:,:,2,:,:,:));%ROn

%imgscan: ACS data that can be used to train GRAPPA coefficients
vrgf_image_raw = twix.image.unsorted();
vrgf_image     = reshape(v'*reshape(vrgf_image_raw,[NFreq_inres Ncoils*NImgLin*NEcho*NSlc*NRep]), sz_image);

%%  the next step is to generate the grappa parameters


vrgf_refscanPC_Part1  = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,1,reordr_P1,:));
vrgf_refscan_Part1    = squeeze(vrgf_refscan(:,:,:,1,reordr_P1,:));
vrgf_phasecor_Part1   = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,1,reordr_P1,:));
vrgf_image_Part1      = squeeze(vrgf_image(:,:,:,1,reordr_P1,:));
refscan_Part1_temp    = PhaseCorrect_yang(vrgf_refscanPC_Part1,vrgf_refscan_Part1);
for i=1:NSeg
  refscan_Part1(:,:,i:NSeg:NRefLin,:)=refscan_Part1_temp(:,:,:,:,i);  
end

sz_k =[NFreq_outres  Ncoils  NImgLin     NSlc   NRep];
Kimage_Part1 = reshape(PhaseCorrect_yang(vrgf_phasecor_Part1,vrgf_image_Part1),sz_k);
disp('---finish phase-correction of Part1--- ');


vrgf_refscanPC_Part2  = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,2,reordr_P2,:));
vrgf_refscan_Part2    = squeeze(vrgf_refscan(:,:,:,2,reordr_P2,:));
vrgf_phasecor_Part2   = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,2,reordr_P2,:));
vrgf_image_Part2      = squeeze(vrgf_image(:,:,:,2,reordr_P2,:));
refscan_Part2_temp    = PhaseCorrect_yang(vrgf_refscanPC_Part2,vrgf_refscan_Part2);
for i=1:NSeg
  refscan_Part2(:,:,i:NSeg:NRefLin,:)=refscan_Part2_temp(:,:,:,:,i);  
end

Kimage_Part2 = reshape(PhaseCorrect_yang(vrgf_phasecor_Part2,vrgf_image_Part2),sz_k);
disp('---finish phase-correction of Part2--- ');

%%  the next step is to generate the grappa kernal
% if ExternalACS is 'true',use the external ACSKernel
if (nargin > 1)
    [~,Ng_external] = recon_std_EPI(varargin{2});
    Ng_Part1 = Ng_external(1:NSlc);
    Ng_Part2 = Ng_external(NSlc+1:2*NSlc);
    
else
    %  the next step is to generate the grappa parameters
    for slc = 1:NSlc
        
        [~,~,Ng_Part1{slc}] = recongrappa_multik_Yang( refscan_Part1(:,:,:,slc), vec(1:NRefLin),'kernel','2x5','dks', R*[1;2] );
        [~,~,Ng_Part2{slc}] = recongrappa_multik_Yang( refscan_Part2(:,:,:,slc), vec(1:NRefLin),'kernel','2x5','dks', R*[1;2] );
    end
end

%% 

% pad the data?
Kimage_Part1_full=zeros(NFreq_outres, Ncoils, NImgLin*R,NSlc,NRep);
Kimage_Part1_full(:,:,R:R:end,:,:) = Kimage_Part1;    



Kimage_Part2_full=zeros(NFreq_outres, Ncoils, NImgLin*R,NSlc,NRep);
Kimage_Part2_full(:,:,R:R:end,:,:) = Kimage_Part2;



%  the next step is to recon
for cntT = 1:NRep
    for slc = 1:NSlc
        
        img_temp1 = recongrappa_multik_Yang( Kimage_Part1_full(:,:,:,slc,cntT),[],'kernel','2x5','N', Ng_Part1{slc}); % important,Ng_...
        img_temp2 = recongrappa_multik_Yang( Kimage_Part2_full(:,:,:,slc,cntT),[],'kernel','2x5','N', Ng_Part2{slc});% important,Ng_...
        
        % correct the partial-Fourier acquiisition, for each coil
        for cntC=1:Ncoils
            img_coil1(:,:,cntC) = flip(reconhd( img_temp1(:,:,cntC),size(img_temp1,1),size(img_temp1,2)));
            img_coil2(:,:,cntC) = flip(reconhd( img_temp2(:,:,cntC),size(img_temp2,1),size(img_temp2,2)));
            
            % flip should be done after reconhd, otherwise it will cause
            % very low SNR for large b-value, e.g. b=2000 or 3000             
            %img_coil1(:,:,cntC) = reconhd( flip(img_temp1(:,:,cntC),1),size(img_temp1,1),size(img_temp1,2) ) ;
            %img_coil2(:,:,cntC) = reconhd( flip(img_temp2(:,:,cntC),1),size(img_temp2,1),size(img_temp2,2) ) ;            
            
        end
        
        % combine the coils to form the final image
        I_Part1(:,:,slc,cntT) = sqrt(sum(abs( fif( img_coil1 ) ).^2,3));
        I_Part2(:,:,slc,cntT) = sqrt(sum(abs( fif( img_coil2 ) ).^2,3));
        
    end
end



[s1, s2, s3, s4]=size(I_Part1);
A=zeros([s1 s2 2*s3 s4]);

A(:,:,1:NSlc,:)        = I_Part1;
A(:,:,NSlc+1:2*NSlc,:) = I_Part2;



