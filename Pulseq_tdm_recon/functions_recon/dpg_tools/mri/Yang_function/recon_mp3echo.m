function [A,B,C] = recon_mp3echo( varargin )
%% add function:SEloop_mode 20200809
twix_  = varargin{1};
twix=twix_{1,2};
hdr = parse_measdat_hdr_Yang(twix.image.filename);
Ncoils = hdr.Ncoils;
NRefLin = hdr.NRefLin;
NSlc = hdr.NSlc;
R = hdr.R;
NEcho = 3;% three echo sequence mp3;
NPhaCorrLin=3;%  3 lines at ky=0 are acquired, RO+ RO- RO+;
NRep = twix.image.NRep;
slic_indexS = hdr.slcindx;
slic_indexM = hdr.slcindx;
slic_indexL = hdr.slcindx;

[~,reordr_S]=sort(slic_indexS);
[~,reordr_M]=sort(slic_indexM);
[~,reordr_L]=sort(slic_indexL);


if (nargin > 1)
   SEloop_mode = varargin{2};   
end

if (nargin > 2)
    twix_acs_ = varargin{3};
    twix_acs  = twix_acs_{1,2};
    hdr_acs = parse_measdat_hdr_Yang(twix_acs.image.filename);
    if ((hdr.R~=hdr_acs.R)||(hdr.NRefLin~=hdr_acs.NRefLin)||(3*hdr.NSlc~=hdr_acs.NSlc))
        error( '---------The dimension of externalACS data does not match------------' );
    end
end


%% % for higher efficiency, high-spatital frequency EPI data is sampled 'on the
% ramps', meaning before the gradient trajectory is stable.  Thus, it
% needs to be 'gridded' onto a constant-velocity-gradient grid prior to
% processing. The VRGF approach does this, e.g. converting 256 sampled-grid
% points into 128 Cartesian-grid points.
%
v = extract_vrgf(twix.image.filename);         % build a VRGF conversion


if (nargin > 2)
    v_acs = extract_vrgf(twix_acs.image.filename);
end

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
refscanPC_raw  = twix.refscanPC.unsorted();
vrgf_refscanPC = reshape(v'*reshape(refscanPC_raw,[NFreq_inres Ncoils*NPhaCorrLin*NEcho*NSlc*NSeg]), sz_refscanPC);
% refscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
vrgf_refscanPC_ROp_ROn(:,:,1,:,:,:) = squeeze((vrgf_refscanPC(:,:,1,:,:,:)+vrgf_refscanPC(:,:,3,:,:,:))/2);%ROp
vrgf_refscanPC_ROp_ROn(:,:,2,:,:,:) = squeeze(vrgf_refscanPC(:,:,2,:,:,:));%ROn


%refscan: ACS data that can be used to train GRAPPA coefficients
vrgf_refscan_raw = twix.refscan.unsorted();
vrgf_refscan = reshape(v'*reshape(vrgf_refscan_raw,[NFreq_inres Ncoils*NRefLin*NEcho*NSlc]), sz_refscan);



%% image data extraction
phasecor_raw   = twix.phasecor.unsorted();
vrgf_phasecor  = reshape(v'*reshape(phasecor_raw,[NFreq_inres Ncoils*NPhaCorrLin*NEcho*NSlc*NRep]), sz_phasecor);
% imgscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','3 echo','NSlc', 'NRep'
vrgf_phasecor_ROp_ROn(:,:,1,:,:,:)= squeeze((vrgf_phasecor(:,:,1,:,:,:)+vrgf_phasecor(:,:,3,:,:,:))/2);%ROp
vrgf_phasecor_ROp_ROn(:,:,2,:,:,:)= squeeze( vrgf_phasecor(:,:,2,:,:,:));%ROn

%imgscan: ACS data that can be used to train GRAPPA coefficients
vrgf_image_raw = twix.image.unsorted();
vrgf_image     = reshape(v'*reshape(vrgf_image_raw,[NFreq_inres Ncoils*NImgLin*NEcho*NSlc*NRep]), sz_image);

%%  the next step is to generate the grappa parameters


%% ########-----data order from mp3echo-----########
%%                            ACSscan     
%       refscanPC:  refscanPC#3   refscanPC#2   refscanPC#1 
%        refscan :     refscan#1     refscan#2     refscan#3
%-->--->--->-->--->--->--sequential data storage-->--->--->-->--->
%
%%                           Imagescan
%        phasecor:    phasecor#3    phasecor#2    phasecor#1 
%           image:       image#1       image#2       image#3
%-->--->--->-->--->--->--sequential data storage-->--->--->-->--->
% The difference between PAT2 and PAT3 could be related to the fact that for PAT3 the PAT reference data are
% usually acquired segmented while for PAT2 this is done in a single shot
%
%% ########---------------------------------########


vrgf_refscanPC_first   = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,3,reordr_S,:));
vrgf_refscan_first     = squeeze(vrgf_refscan(:,:,:,1,reordr_S,:));
vrgf_phasecor_first    = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,3,reordr_S,:));
vrgf_image_first       = squeeze(vrgf_image(:,:,:,1,reordr_S,:));
refscan_first_temp     = PhaseCorrect_yang(vrgf_refscanPC_first,vrgf_refscan_first);
for i=1:NSeg
  refscan_first(:,:,i:NSeg:NRefLin,:)=refscan_first_temp(:,:,:,:,i);  
end

sz_k =[NFreq_outres  Ncoils  NImgLin     NSlc   NRep];
Kimage_first = reshape(PhaseCorrect_yang(vrgf_phasecor_first,vrgf_image_first),sz_k);
disp('---finish phase-correction of first  TE--- ');


vrgf_refscanPC_second  = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,2,reordr_M,:));
vrgf_refscan_second    = squeeze(vrgf_refscan(:,:,:,2,reordr_M,:));
vrgf_phasecor_second   = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,2,reordr_M,:));
vrgf_image_second      = squeeze(vrgf_image(:,:,:,2,reordr_M,:));
refscan_second_temp    = PhaseCorrect_yang(vrgf_refscanPC_second,vrgf_refscan_second);
for i=1:NSeg
  refscan_second(:,:,i:NSeg:NRefLin,:)=refscan_second_temp(:,:,:,:,i);  
end

Kimage_second = reshape(PhaseCorrect_yang(vrgf_phasecor_second,vrgf_image_second),sz_k);
disp('---finish phase-correction of second TE--- ');



vrgf_refscanPC_third   = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,1,reordr_L,:));
vrgf_refscan_third     = squeeze(vrgf_refscan(:,:,:,3,reordr_L,:));
vrgf_phasecor_third    = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,1,reordr_L,:));
vrgf_image_third       = squeeze(vrgf_image(:,:,:,3,reordr_L,:));
refscan_third_temp     = PhaseCorrect_yang(vrgf_refscanPC_third,vrgf_refscan_third);
for i=1:NSeg
  refscan_third(:,:,i:NSeg:NRefLin,:)=refscan_third_temp(:,:,:,:,i);  
end
 
Kimage_third = reshape(PhaseCorrect_yang(vrgf_phasecor_third,vrgf_image_third),sz_k);
disp('---finish phase-correction of third  TE--- ');




%%  the next step is to generate the grappa kernel
% if ExternalACS is 'true',use the external ACSKernel
if (nargin > 2)
    [~,Ng_external] = recon_std_EPI(varargin{3});
%     Ng_third = Ng_external;
%     Ng_first = Ng_external(NSlc+1:2*NSlc);
%     Ng_second = Ng_external(2*NSlc+1:3*NSlc);

    k1 = Ng_external;
    k2 = Ng_external(NSlc+1:2*NSlc);
    k3 = Ng_external(2*NSlc+1:3*NSlc);
    if SEloop_mode==0
        Ng_third = k1;
        Ng_first = k2;
        Ng_second = k3;
        
    elseif SEloop_mode==1
        Ng_second = k1;
        Ng_third = k2;
        Ng_first = k3;
    elseif SEloop_mode==2
        Ng_first = k1;
        Ng_second = k2;
        Ng_third = k3;    
    end
      
else
    
    % use the internal ACSKernel
    for slc = 1:NSlc
        [~,~,Ng_first{slc}]  = recongrappa_multik_Yang( refscan_first(:,:,:,slc),  vec(1:NRefLin),'kernel','2x5','dks', R*[1;2] );
        [~,~,Ng_second{slc}] = recongrappa_multik_Yang( refscan_second(:,:,:,slc), vec(1:NRefLin),'kernel','2x5','dks', R*[1;2] );
        [~,~,Ng_third{slc}]  = recongrappa_multik_Yang( refscan_third(:,:,:,slc),  vec(1:NRefLin),'kernel','2x5','dks', R*[1;2] );
    end
       
end

 

%% 

% pad the data?
Kimage_first_full=zeros(NFreq_outres, Ncoils, NImgLin*R,NSlc,NRep);
Kimage_first_full(:,:,R:R:end,:,:) = Kimage_first;    

Kimage_second_full=zeros(NFreq_outres, Ncoils, NImgLin*R,NSlc,NRep);
Kimage_second_full(:,:,R:R:end,:,:) = Kimage_second;

Kimage_third_full=zeros(NFreq_outres, Ncoils, NImgLin*R,NSlc,NRep);
Kimage_third_full(:,:,R:R:end,:,:) = Kimage_third;




%  the next step is to recon
for cntT = 1:NRep
% for cntT = 1:1
    for slc = 1:NSlc
        
        img_temp1 = recongrappa_multik_Yang( Kimage_first_full(:,:,:,slc,cntT), [],'kernel','2x5','N',Ng_first{slc}); % important,Ng_...
        img_temp2 = recongrappa_multik_Yang( Kimage_second_full(:,:,:,slc,cntT),[],'kernel','2x5','N',Ng_second{slc});% important,Ng_...
        img_temp3 = recongrappa_multik_Yang( Kimage_third_full(:,:,:,slc,cntT), [],'kernel','2x5','N',Ng_third{slc});% important,Ng_...
        % correct the partial-Fourier acquiisition, for each coil
        for cntC=1:Ncoils
            img_coil1(:,:,cntC) = flip(reconhd( img_temp1(:,:,cntC),size(img_temp1,1),size(img_temp1,2)));
            img_coil2(:,:,cntC) = flip(reconhd( img_temp2(:,:,cntC),size(img_temp2,1),size(img_temp2,2)));
            img_coil3(:,:,cntC) = flip(reconhd( img_temp3(:,:,cntC),size(img_temp3,1),size(img_temp3,2)));

            % flip should be done after reconhd, otherwise it will cause
            % very low SNR for large b-value, e.g. b=2000 or 3000              
            %img_coil1(:,:,cntC) = reconhd( flip(img_temp1(:,:,cntC),1),size(img_temp1,1),size(img_temp1,2) ) ;
            %img_coil2(:,:,cntC) = reconhd( flip(img_temp2(:,:,cntC),1),size(img_temp2,1),size(img_temp2,2) ) ;
            %img_coil3(:,:,cntC) = reconhd( flip(img_temp3(:,:,cntC),1),size(img_temp3,1),size(img_temp3,2) ) ;            
            
        end
        
        % combine the coils to form the final image
        I_short(:,:,slc,cntT)  = sqrt(sum(abs( fif( img_coil1 ) ).^2,3));
        I_medium(:,:,slc,cntT) = sqrt(sum(abs( fif( img_coil2 ) ).^2,3));
        I_long(:,:,slc,cntT)   = sqrt(sum(abs( fif( img_coil3 ) ).^2,3));
        cntT
        slc
    end
	
end

A=I_short;
B=I_medium;
C=I_long;

