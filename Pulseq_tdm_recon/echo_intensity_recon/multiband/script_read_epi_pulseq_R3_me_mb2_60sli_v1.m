close all; clear all ; clc
addpath(genpath('../../functions_recon'))
% addpath(genpath('/data/pnl/home/ql087/freesurfer'))
addpath(genpath('../../Pulseq_Qiang_v142-main'))
% addpath(genpath('/data/pnl/home/ql087/Bruker_2022'))
% addpath(genpath('/data/pnl/home/ql087/VecNorm'))

save_path='/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/data_processing/TDM_revision/sub_3/';
save_filename = 'pulseq_me_te456.mat';

data_path = '/data/pnlx/home/ln915/Data/TDMvsME/sub3/raw/'
filename1 = 'meas_MID01198_FID61553_pulseq_mb2_me_te456_a';
filename2 = 'meas_MID01199_FID61554_pulseq_mb2_me_te456_b';
filename3 = 'meas_MID01200_FID61555_pulseq_mb2_me_te456_c';
filename4 = 'meas_MID01201_FID61556_pulseq_60sli_3_shot_ref'; % ACS

%% section a
D=dir([data_path filename1]);
[~,I]=sort([D(:).datenum]);
twix_obj = mapVBVD([data_path filename1]);
seq = mr.Sequence();
read(seq,'/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Tests/Test_2023/Test_72_Dec_26/signleecho_mb1_30sli_2p5mm_R3_multib_te1_recon.seq')
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%% define raw data

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end
rawdata=single(rawdata);
rawdata=reshape(rawdata,[212 32 22 3 10 28 3]);
rawdata=permute(rawdata,[1 2 3 5 6 7 4]);
rawdata=reshape(rawdata,[212 32 55440]);


traj_recon_delay=0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);

%% manually (re-)define FOV and resolution, I found we need it for no blip EPI, Qiang Liu
Nx=88; Ny=Nx; Ny_sampled=ceil(Ny*6/8/3);
%% analyze the trajecotory, resample the data
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions

ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

slice_num=30;
Ngroup=slice_num;
slice_num_per_group=10;
dire=28;
bvalue=3;
diff_plus_dum=dire*bvalue; 
SMS=2;
evp.NSlcMeas=SMS*slice_num;
multiply=slice_num_per_group*diff_plus_dum;
multipl_dum=slice_num_per_group*1; %3 TEs ?
nAd=nAcq/diff_plus_dum;
data_resampled=zeros(length(kxx),nCoils, nAd, diff_plus_dum);
% ktraj_adc2=reshape(ktraj_adc2,[size(ktraj_adc2,1) size(ktraj_adc2,2) size(ktraj_adc2,3)/3/5 3*5]);
rawdata=reshape(rawdata,[size(rawdata,1)  size(rawdata,2) nAd diff_plus_dum]);

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

tic
parfor t_c=1:size(data_resampled,4)
    for a=1:nAd
        for c=1:nCoils
            data_resampled(:,c,a,t_c)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a, t_c),kxx,'spline',0);
        end
    end
end
toc

delete(gcp('nocreate'))

data_resampled=reshape(data_resampled,[length(kxx), nCoils, nAcq]);

%% save the reference data
ref=permute(data_resampled,[1 3 2]); % 64, 32, 168
ref=reshape(ref,[Nx Ny_sampled multiply*3 nCoils]);

%% start to do LPC for EPI
diff_plus_dum=diff_plus_dum*3;
k=ref(:,:,(multipl_dum+1):end,:);
k=reshape(k,[Nx Ny_sampled multipl_dum diff_plus_dum-1 nCoils]);
k_gc=zeros(Ny_sampled,Nx,multipl_dum,diff_plus_dum-1,nCoils);
for dti_loop=1:(diff_plus_dum-1)

    ref=ref(:,:,1:multipl_dum,:);

    for i=1:size(ref,3)

        s=squeeze(ref(:,8:10,i,:));

        k_1=squeeze(k(:,:,i,dti_loop,:));
        %% LPC based on WSH's LPC code
        S0 = fif(mean(s(:,[1 3],:),2));
        S1 = fif(mean(s(:,[ 2 ],:),2));
        for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
        k_gc(:,:,i,dti_loop,:) = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);

    end
end

k_gc=squeeze(reshape(k_gc,[Ny_sampled, Nx, slice_num_per_group, (diff_plus_dum-1), nCoils]));
k_gc=permute(k_gc,[2 5 1 3 4]);% 110 32 28 2556
sz_k=[Nx nCoils Ny_sampled slice_num_per_group (diff_plus_dum-1)]; %I treat the gSlider loop as the diffusion
Kimage_short = reshape(k_gc,sz_k);
Kimage_short_full=zeros(Nx,nCoils, Ny_sampled, slice_num, (diff_plus_dum-1));
Kimage_short_full(:,:,:,1:3:end,:)=Kimage_short;
clear rawdata data_resampled k_gc Kimage_short

%% section b
D=dir([data_path filename2]);
[~,I]=sort([D(:).datenum]);
twix_obj = mapVBVD([data_path filename2]);
seq = mr.Sequence();
read(seq,'/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Tests/Test_2023/Test_72_Dec_26/signleecho_mb1_30sli_2p5mm_R3_multib_te1_recon.seq')
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%% define raw data

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end
rawdata=single(rawdata);
rawdata=reshape(rawdata,[212 32 22 3 10 28 3]);
rawdata=permute(rawdata,[1 2 3 5 6 7 4]);
rawdata=reshape(rawdata,[212 32 55440]);


traj_recon_delay=0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);

%% manually (re-)define FOV and resolution, I found we need it for no blip EPI, Qiang Liu
Nx=88; Ny=Nx; Ny_sampled=ceil(Ny*6/8/3);
%% analyze the trajecotory, resample the data
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions

ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

slice_num=30;
Ngroup=slice_num;
slice_num_per_group=10;
dire=28;
bvalue=3;
diff_plus_dum=dire*bvalue; 
SMS=2;
evp.NSlcMeas=SMS*slice_num;
multiply=slice_num_per_group*diff_plus_dum;
multipl_dum=slice_num_per_group*1; %3 TEs ?
nAd=nAcq/diff_plus_dum;
data_resampled=zeros(length(kxx),nCoils, nAd, diff_plus_dum);
% ktraj_adc2=reshape(ktraj_adc2,[size(ktraj_adc2,1) size(ktraj_adc2,2) size(ktraj_adc2,3)/3/5 3*5]);
rawdata=reshape(rawdata,[size(rawdata,1)  size(rawdata,2) nAd diff_plus_dum]);

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

tic
parfor t_c=1:size(data_resampled,4)
    for a=1:nAd
        for c=1:nCoils
            data_resampled(:,c,a,t_c)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a, t_c),kxx,'spline',0);
        end
    end
end
toc

delete(gcp('nocreate'))

data_resampled=reshape(data_resampled,[length(kxx), nCoils, nAcq]);

%% save the reference data
ref=permute(data_resampled,[1 3 2]); % 64, 32, 168
ref=reshape(ref,[Nx Ny_sampled multiply*3 nCoils]);

%% start to do LPC for EPI
diff_plus_dum=diff_plus_dum*3;
k=ref(:,:,(multipl_dum+1):end,:);
k=reshape(k,[Nx Ny_sampled multipl_dum diff_plus_dum-1 nCoils]);
k_gc=zeros(Ny_sampled,Nx,multipl_dum,diff_plus_dum-1,nCoils);
for dti_loop=1:(diff_plus_dum-1)

    ref=ref(:,:,1:multipl_dum,:);

    for i=1:size(ref,3)

        s=squeeze(ref(:,8:10,i,:));

        k_1=squeeze(k(:,:,i,dti_loop,:));
        %% LPC based on WSH's LPC code
        S0 = fif(mean(s(:,[1 3],:),2));
        S1 = fif(mean(s(:,[ 2 ],:),2));
        for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
        k_gc(:,:,i,dti_loop,:) = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);

    end
end

k_gc=squeeze(reshape(k_gc,[Ny_sampled, Nx, slice_num_per_group, (diff_plus_dum-1), nCoils]));
k_gc=permute(k_gc,[2 5 1 3 4]);% 110 32 28 2556
sz_k=[Nx nCoils Ny_sampled slice_num_per_group (diff_plus_dum-1)]; %I treat the gSlider loop as the diffusion
Kimage_short = reshape(k_gc,sz_k);
Kimage_short_full(:,:,:,2:3:end,:)=Kimage_short;
clear rawdata data_resampled k_gc Kimage_short

%% section c
D=dir([data_path filename3]);
[~,I]=sort([D(:).datenum]);
twix_obj = mapVBVD([data_path filename3]);
seq = mr.Sequence();
read(seq,'/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Tests/Test_2023/Test_72_Dec_26/signleecho_mb1_30sli_2p5mm_R3_multib_te1_recon.seq')
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%% define raw data

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end
rawdata=single(rawdata);
rawdata=reshape(rawdata,[212 32 22 3 10 28 3]);
rawdata=permute(rawdata,[1 2 3 5 6 7 4]);
rawdata=reshape(rawdata,[212 32 55440]);


traj_recon_delay=0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);

%% manually (re-)define FOV and resolution, I found we need it for no blip EPI, Qiang Liu
Nx=88; Ny=Nx; Ny_sampled=ceil(Ny*6/8/3);
%% analyze the trajecotory, resample the data
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions

ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

slice_num=30;
Ngroup=slice_num;
slice_num_per_group=10;
dire=28;
bvalue=3;
diff_plus_dum=dire*bvalue; 
SMS=2;
evp.NSlcMeas=SMS*slice_num;
multiply=slice_num_per_group*diff_plus_dum;
multipl_dum=slice_num_per_group*1; %3 TEs ?
nAd=nAcq/diff_plus_dum;
data_resampled=zeros(length(kxx),nCoils, nAd, diff_plus_dum);
% ktraj_adc2=reshape(ktraj_adc2,[size(ktraj_adc2,1) size(ktraj_adc2,2) size(ktraj_adc2,3)/3/5 3*5]);
rawdata=reshape(rawdata,[size(rawdata,1)  size(rawdata,2) nAd diff_plus_dum]);

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

tic
parfor t_c=1:size(data_resampled,4)
    for a=1:nAd
        for c=1:nCoils
            data_resampled(:,c,a,t_c)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a, t_c),kxx,'spline',0);
        end
    end
end
toc

delete(gcp('nocreate'))

data_resampled=reshape(data_resampled,[length(kxx), nCoils, nAcq]);

%% save the reference data
ref=permute(data_resampled,[1 3 2]); % 64, 32, 168
ref=reshape(ref,[Nx Ny_sampled multiply*3 nCoils]);

%% start to do LPC for EPI
diff_plus_dum=diff_plus_dum*3;
k=ref(:,:,(multipl_dum+1):end,:);
k=reshape(k,[Nx Ny_sampled multipl_dum diff_plus_dum-1 nCoils]);
k_gc=zeros(Ny_sampled,Nx,multipl_dum,diff_plus_dum-1,nCoils);
for dti_loop=1:(diff_plus_dum-1)

    ref=ref(:,:,1:multipl_dum,:);

    for i=1:size(ref,3)

        s=squeeze(ref(:,8:10,i,:));

        k_1=squeeze(k(:,:,i,dti_loop,:));
        %% LPC based on WSH's LPC code
        S0 = fif(mean(s(:,[1 3],:),2));
        S1 = fif(mean(s(:,[ 2 ],:),2));
        for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
        k_gc(:,:,i,dti_loop,:) = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);

    end
end

k_gc=squeeze(reshape(k_gc,[Ny_sampled, Nx, slice_num_per_group, (diff_plus_dum-1), nCoils]));
k_gc=permute(k_gc,[2 5 1 3 4]);% 110 32 28 2556
sz_k=[Nx nCoils Ny_sampled slice_num_per_group (diff_plus_dum-1)]; %I treat the gSlider loop as the diffusion
Kimage_short = reshape(k_gc,sz_k);
Kimage_short_full(:,:,:,3:3:end,:)=Kimage_short;
clear rawdata data_resampled k_gc Kimage_short



%% ACS data
D=dir([data_path filename4]);
[~,I]=sort([D(:).datenum]);
twix_obj = mapVBVD([data_path filename4]);
seq = mr.Sequence();
read(seq,'/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Tests/Test_2023/Test_72_Dec_26/epidiff_3_shot_ref_2p5mm_60sli_appa_recon.seq') % I want to keep this line, QL
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end
traj_recon_delay=0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);

%% manually (re-)define FOV and resolution, I found we need it for no blip EPI, Qiang Liu
Nx=88; Ny=Nx; Ny_sampled=ceil(Ny*6/8/3);
%% analyze the trajecotory, resample the data
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);
slice_num=SMS*slice_num;
kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions
ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

data_resampled=zeros(length(kxx), nCoils, nAcq);
for a=1:nAcq
    for c=1:nCoils
        data_resampled(:,c,a)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a),kxx,'spline',0);
    end
end

ref=permute(data_resampled,[1 3 2]); % 64, 32, 168
ref=ref(:,1:size(ref,2)*2/3,:);
diff_plus_dum=3;
ref=reshape(ref,[Nx Ny_sampled slice_num*3*(diff_plus_dum-1) nCoils]);
multiply_ref=slice_num*diff_plus_dum;
R=3;
%% start to do LPC for EPI
k=ref(:,:,(multiply_ref+1):end,:);
k_gc_1=zeros(Ny_sampled,Nx,multiply_ref,nCoils);
ref=ref(:,:,1:multiply_ref,:);

for i=1:size(ref,3)

    s=squeeze(ref(:,2:4,i,:));
    k_1=squeeze(k(:,:,i,:));
    %% LPC based on WSH's LPC code
    S0 = fif(mean(s(:,[1 3],:),2));
    S1 = fif(mean(s(:,[ 2 ],:),2));
    for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
    k_gc_1(:,:,i,:) = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);

end

clear rawdata data_resampled
load('/data/pnl/home/ql087/data_bwh/prot_mb4.mat')
PhaseShift=2*pi/2;


acs=permute(k_gc_1(2:13,:,:,:), [2 4 1 3]);
acs=reshape(acs,[Nx nCoils 12 slice_num R]);
refscan_first=zeros(Nx,nCoils,12*R,slice_num);
for iii=1:R
    refscan_first(:,:,iii:R:end,:)=squeeze(acs(:,:,:,:,iii));
end



k_trgt=permute(refscan_first,[1 3 4 2]);
p=[2:2:slice_num,1:2:slice_num];
[~,porder]=sort(p);
k_trgt=k_trgt(:,:,porder,:);
% k_trgt=k_trgt(:,:,end:-1:1,:);
clear refscan_first acs

sz = size(k_trgt);

if ~exist('slices','var'); slices = 1:slice_num; end;
if ~exist('chi','var'); chi = 1e-6; end;
if ~exist('eta','var'); eta = 1; end;


phzabs = exp(j*(PhaseShift*(0:(size(k_trgt,2)-1))/R - (0.0*pi) - PhaseShift ) );

for cnt=0:(SMS-1),
    k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
end;



if ~exist('Npppp','var')

    fprintf('%d',mod(1:length(slices),10));fprintf('\n');
    for slc = 1:slice_num;
        fprintf('o');
        [~,~,Np{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',R*[1 2],...
            'chi',chi,'eta',eta );
    end;
    fprintf('\n');

end

z = [];
for cnt=1:length(Np{slices(1)})
    dk = diff(find( Np{slices(1)}(cnt).pattern == '*' ));
    z(dk) =  1;
end;

%% recon

%%%%calculate w, repeat Ngroup time
for slc = 1:Ngroup
    iii = 1:size(k_trgt,2);
    sz = size(k_trgt);
    for cnt=1:SMS,
        in2{cnt} = zeros( sz([1 2 4]) );
        in2{cnt}(:,(sz(2)-length(iii))/2+iii,:) = squeeze( k_trgt(:,:,slc+(cnt-1)*Ngroup,:) );
    end;
    [~,w{slc}] =  MultisliceGRAPPA_2kernal_leakBlock( in2{1}, in2, [ 5 3 1 R ], 'full', prot);

end
%--------------------------------------------------------------------

NRep=size(Kimage_short_full,5);
k_data_gc0 = permute(Kimage_short_full,[1 3 4 2 5]);


for cRep=1:NRep

    %     k_data_gc0 = permute(Kimage_short,[1 3 4 2 5]);
    k_data_gc  = k_data_gc0(:,:,:,:,cRep);
    k_data_gc_deblur = CaipirinhaDeblur_v3_wsh_ql( k_data_gc(:,:,:,:) );

    if exist('runmod')
        if ( runmod == 1 )
            k_data_gc_deblur = k_data_gc_deblur(:,:,:,coils);
        end;
    end;

    nky = size(k_data_gc_deblur,2);


    fprintf('%d',mod(1:length(slices),10)); fprintf('\n');
    for slc = 1:Ngroup;
        fprintf('.');

        in1 = zeros([ size(k_data_gc_deblur,1) nky*R sz(4)]);

        in1(:,R*(0:size(k_data_gc_deblur,2)-1)+1,:) = squeeze( k_data_gc_deblur(:,:,slc,:) );

        tmp = squeeze(MultisliceGRAPPA_2kernal_leakBlock( in1, w{slc}, [ 5 3 1 R ]) );


        phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/R - (0.0*pi) - PhaseShift ) );
        slcgrp = slc + [ 0:(SMS-1) ]*Ngroup;

        sz_in1 = size(in1);
        for cnt=1:SMS,
            curslc = slcgrp(cnt);
            Fa1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,1:2*R:end,:,cnt),[2 1 3]),1:2*R:sz_in1(2),'kernel','2x5','N',Np{curslc});
            Fb1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,(1+R):2*R:end,:,cnt),[2 1 3]),(1+R):2*R:sz_in1(2),'kernel','2x5','N',Np{curslc});

            Fa2 = permute(tunfold(fif(Fa1),2),[2 1]);
            Fb2 = permute(tunfold(fif(Fb1),2),[2 1]);

            Fb3 = phzshift( Fa2, Fb2,{'nofft','nocombo'} );
            Fb4 = ifi( trefold(permute(Fb3,[2 1]),sz_in1([2 1 3]),2) );

            Fc1 = zeros(size(Fa1));
            Fc1( 1:2*R:end, :, : ) = Fa1( 1:2*R:end, :, : );
            Fc1( (1+R):2*R:end, :, : ) = Fb4( (1+R):2*R:end, :, : );

            Fd1 = recongrappa_multik(size(Fc1),Fc1,[],'kernel','2x5','dks',R,'N',Np{curslc});

            F2klb(:,:,curslc,:) = Fd1;
            F2klb(:,:,curslc,:) =  tmult( F2klb(:,:,curslc,:), diag(conj(phzabs).^(cnt-1)), 1);
        end;

        %if verbose, keyboard; end;
    end
    fprintf('\n');
    k_all(:,:,:,:,cRep)=F2klb;

end

%% partial Fourier+

aa=k_all;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Pulseq_Mprage_Recon_Toolbox'))

for cntsli=1:size(aa,3)

    for cntcoil=1:size(aa,4)
        for cRep=1:size(aa,5)
            img_aa(:,:,cntcoil,cRep)=flip(reconhd(aa(:,:,cntsli,cntcoil,cRep),size(aa,1),size(aa,2)));
            img_aa(:,:,cntcoil,cRep)=apodize3d_QL_v2(img_aa(:,:,cntcoil, cRep), 0.25);
        end

    end
    img(:,:,cntsli,:)=sqrt(sum(abs( fif( img_aa ) ).^2,3));
end
img_final=zeros(88,88,60,252);
img_final(:,:,:,2:end)=img;
clear img
img_final=reshape(img_final,[88,88,60,28,3,3]);
tmp=[2 4:28];
img_final=img_final(:,:,:,tmp,:,:);
clear tmp

save([save_path, save_filename],'img_final','-v7.3') 
