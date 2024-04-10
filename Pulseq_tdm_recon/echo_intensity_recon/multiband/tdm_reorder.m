clear;close all;clc;
% note that tdm 456/123

Base_dir='/data/pnl/home/ql087/data_processing/2024_01_05_pulseq_multiecho_mb_invivo/invivo_tdm/processing/';

%% 
tmp1=strcat(Base_dir, 'pulseq_tdm_4.mat');
tmp2=strcat(Base_dir, 'pulseq_tdm_5.mat');
tmp3=strcat(Base_dir, 'pulseq_tdm_6.mat');
load(tmp1);
tdm1=img_final;
load(tmp2)
tdm2=img_final;
load(tmp3)
tdm3=img_final;
clear img_final tmp*

b_num=3;
dir=26;
tmp1=tdm1;
tmp2=tdm2;
tmp3=tdm3;
tdm1=zeros(size(tmp1));
tdm2=tdm1;tdm3=tdm1;

tdm1(:,:,1:10,:,:)=tmp1(:,:,1:10,:,:);
tdm1(:,:,11:10*2,:,:)=tmp2(:,:,11:10*2,:,:);
tdm1(:,:,21:10*3,:,:)=tmp3(:,:,21:10*3,:,:);

tdm1(:,:,31:40,:,:)=tmp1(:,:,31:40,:,:);
tdm1(:,:,41:10*5,:,:)=tmp2(:,:,41:10*5,:,:);
tdm1(:,:,51:10*6,:,:)=tmp3(:,:,51:10*6,:,:);


tdm2(:,:,1:10,:,:)=tmp3(:,:,1:10,:,:);
tdm2(:,:,11:10*2,:,:)=tmp1(:,:,11:10*2,:,:);
tdm2(:,:,21:10*3,:,:)=tmp2(:,:,21:10*3,:,:);

tdm2(:,:,31:40,:,:)=tmp3(:,:,31:40,:,:);
tdm2(:,:,41:10*5,:,:)=tmp1(:,:,41:10*5,:,:);
tdm2(:,:,51:10*6,:,:)=tmp2(:,:,51:10*6,:,:);

tdm3(:,:,1:10,:,:)=tmp2(:,:,1:10,:,:);
tdm3(:,:,11:10*2,:,:)=tmp3(:,:,11:10*2,:,:);
tdm3(:,:,21:10*3,:,:)=tmp1(:,:,21:10*3,:,:);

tdm3(:,:,31:40,:,:)=tmp2(:,:,31:40,:,:);
tdm3(:,:,41:10*5,:,:)=tmp3(:,:,41:10*5,:,:);
tdm3(:,:,51:10*6,:,:)=tmp1(:,:,51:10*6,:,:);

clear tmp*
tdm4=tdm1;
tdm5=tdm2;
tdm6=tdm3;
save('tdm_te456.mat','tdm4','tdm5','tdm6')
% save('tdm_te123.mat','tdm1','tdm2','tdm3')