%% mp3echo Part GRAPPA R =  3
clc;clear;close all;
directory = '/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/data_bwh_mgh/2023_04_03_bwh_pulseq_phantom/IN_VIVO/';  

%% 

s1 = 'meas_MID00066_FID857184_mb_mp3e_L0';
s2 = 'meas_MID00072_FID857190_mb_mp3e_L1';
s3 = 'meas_MID00076_FID857194_mb_mp3e_L2';
s4 = 'meas_MID00059_FID857177_ACS_STD_mp3e_slice60';

k_mp_L0  = mapVBVD(   findfile(directory,s1)   );
k_mp_L1  = mapVBVD(   findfile(directory,s2)   );
k_mp_L2  = mapVBVD(   findfile(directory,s3)   );
k_std    = mapVBVD(   findfile(directory,s4)   );

[P2_TE1,P3_TE2,P1_TE3] = recon_sms_mp3echo(k_mp_L0,0,k_std);
[P3_TE1,P1_TE2,P2_TE3] = recon_sms_mp3echo(k_mp_L1,1,k_std);
save('L1.mat','P3_TE1', 'P1_TE2', 'P2_TE3')
[P1_TE1,P2_TE2,P3_TE3] = recon_sms_mp3echo(k_mp_L2,2,k_std);
save('L2.mat','P1_TE1', 'P2_TE2', 'P3_TE3')

img1_TE1 = cat(3,P1_TE1,P2_TE1,P3_TE1);
img1_TE2 = cat(3,P1_TE2,P2_TE2,P3_TE2);
img1_TE3 = cat(3,P1_TE3,P2_TE3,P3_TE3);

save([directory 'mb_mp.mat']);
%% 

% std_mb1  = 'mb_std_Part1_TE79';
% std_mb2  = 'mb_std_Part2_TE79';
% std_mb3  = 'mb_std_Part3_TE79';
% std_mb4  = 'mb_std_Part1_TE108';
% std_mb5  = 'mb_std_Part2_TE108';
% std_mb6  = 'mb_std_Part3_TE108';
% std_mb7  = 'mb_std_Part1_TE137';
% std_mb8  = 'mb_std_Part2_TE137';
% std_mb9  = 'mb_std_Part3_TE137';
% k_stdMB_Part1_TE1 = mapVBVD(   findfile(directory,std_mb1)   );
% k_stdMB_Part2_TE1 = mapVBVD(   findfile(directory,std_mb2)   );
% k_stdMB_Part3_TE1 = mapVBVD(   findfile(directory,std_mb3)   );
% 
% k_stdMB_Part1_TE2 = mapVBVD(   findfile(directory,std_mb4)   );
% k_stdMB_Part2_TE2 = mapVBVD(   findfile(directory,std_mb5)   );
% k_stdMB_Part3_TE2 = mapVBVD(   findfile(directory,std_mb6)  );
% 
% k_stdMB_Part1_TE3 = mapVBVD(   findfile(directory,std_mb7)  );
% k_stdMB_Part2_TE3 = mapVBVD(   findfile(directory,std_mb8)  );
% k_stdMB_Part3_TE3 = mapVBVD(   findfile(directory,std_mb9)  );
% 
% [stdMB_P1_TE1]  = recon_sms_std(k_stdMB_Part1_TE1);
% [stdMB_P2_TE1]  = recon_sms_std(k_stdMB_Part2_TE1);
% [stdMB_P3_TE1]  = recon_sms_std(k_stdMB_Part3_TE1);
% 
% [stdMB_P1_TE2]  = recon_sms_std(k_stdMB_Part1_TE2);
% [stdMB_P2_TE2]  = recon_sms_std(k_stdMB_Part2_TE2);
% [stdMB_P3_TE2]  = recon_sms_std(k_stdMB_Part3_TE2);
% 
% [stdMB_P1_TE3]  = recon_sms_std(k_stdMB_Part1_TE3);
% [stdMB_P2_TE3]  = recon_sms_std(k_stdMB_Part2_TE3);
% [stdMB_P3_TE3]  = recon_sms_std(k_stdMB_Part3_TE3);
% 
% 
% img2_TE1 = cat(3,stdMB_P1_TE1,stdMB_P2_TE1,stdMB_P3_TE1);
% img2_TE2 = cat(3,stdMB_P1_TE2,stdMB_P2_TE2,stdMB_P3_TE2);
% img2_TE3 = cat(3,stdMB_P1_TE3,stdMB_P2_TE3,stdMB_P3_TE3);
% 
% save([directory 'mb.mat']);
%% 

% std_sb1  = 'sb_std_Part1_TE79';
% std_sb2  = 'sb_std_Part2_TE79';
% std_sb3  = 'sb_std_Part3_TE79';
% std_sb4  = 'sb_std_Part1_TE108';
% std_sb5  = 'sb_std_Part2_TE108';
% std_sb6  = 'sb_std_Part3_TE108';
% std_sb7  = 'sb_std_Part1_TE137';
% std_sb8  = 'sb_std_Part2_TE137';
% std_sb9  = 'sb_std_Part3_TE137';
% k_stdSB_Part1_TE1 = mapVBVD(   findfile(directory,std_sb1)   );
% k_stdSB_Part2_TE1 = mapVBVD(   findfile(directory,std_sb2)   );
% k_stdSB_Part3_TE1 = mapVBVD(   findfile(directory,std_sb3)   );
% 
% k_stdSB_Part1_TE2 = mapVBVD(   findfile(directory,std_sb4)   );
% k_stdSB_Part2_TE2 = mapVBVD(   findfile(directory,std_sb5)   );
% k_stdSB_Part3_TE2 = mapVBVD(   findfile(directory,std_sb6)  );
% 
% k_stdSB_Part1_TE3 = mapVBVD(   findfile(directory,std_sb7)  );
% k_stdSB_Part2_TE3 = mapVBVD(   findfile(directory,std_sb8)  );
% k_stdSB_Part3_TE3 = mapVBVD(   findfile(directory,std_sb9)  );
% 
% [stdSB_P1_TE1,~]  = recon_std_EPI(k_stdSB_Part1_TE1);
% [stdSB_P2_TE1,~]  = recon_std_EPI(k_stdSB_Part2_TE1);
% [stdSB_P3_TE1,~]  = recon_std_EPI(k_stdSB_Part3_TE1);
% 
% [stdSB_P1_TE2,~]  = recon_std_EPI(k_stdSB_Part1_TE2);
% [stdSB_P2_TE2,~]  = recon_std_EPI(k_stdSB_Part2_TE2);
% [stdSB_P3_TE2,~]  = recon_std_EPI(k_stdSB_Part3_TE2);
% 
% [stdSB_P1_TE3,~]  = recon_std_EPI(k_stdSB_Part1_TE3);
% [stdSB_P2_TE3,~]  = recon_std_EPI(k_stdSB_Part2_TE3);
% [stdSB_P3_TE3,~]  = recon_std_EPI(k_stdSB_Part3_TE3);
% 
% img3_TE1 = cat(3,stdSB_P1_TE1,stdSB_P2_TE1,stdSB_P3_TE1);
% img3_TE2 = cat(3,stdSB_P1_TE2,stdSB_P2_TE2,stdSB_P3_TE2);
% img3_TE3 = cat(3,stdSB_P1_TE3,stdSB_P2_TE3,stdSB_P3_TE3);
% 
% save([directory 'sb.mat']);
%% 
% 
% ind=4;
% figure;
% subplot(311);imshow3(   mp_TE1(:,:,[21:40],ind),[ ],[2 10]);
% title('TDM-SMS ', 'FontSize', 10);
% subplot(312);imshow3(stdMB_P2_TE1(:,:,:,ind),[ ],[2 10]);
% title('MB ', 'FontSize', 10);
% subplot(313);imshow3(stdSB_P2_TE1(:,:,:,ind),[ ],[2 10]);
% title('SB ', 'FontSize', 10);
% 
% figure;
% subplot(311);imshow3(   mp_TE2(:,:,[21:40],ind),[],[2 10]);
% title('TDM-SMS ', 'FontSize', 10);
% subplot(312);imshow3(stdMB_P2_TE2(:,:,:,ind),[],[2 10]);
% title('MB ', 'FontSize', 10);
% subplot(313);imshow3(stdSB_P2_TE2(:,:,:,ind),[],[2 10]);
% title('SB ', 'FontSize', 10);
% 
% figure;
% subplot(311);imshow3(   mp_TE3(:,:,[21:40],ind),[],[2 10]);
% title('TDM-SMS ', 'FontSize', 10);
% subplot(312);imshow3(stdMB_P2_TE3(:,:,:,ind),[],[2 10]);
% title('MB ', 'FontSize', 10);
% subplot(313);imshow3(stdSB_P2_TE3(:,:,:,ind),[],[2 10]);
% title('SB ', 'FontSize', 10);
% 
% 
% 
% % a=2e-05/4;
% % figure;imshow3(mp_long,[0 a],[6 10]);
% % figure;imshow3(std_TE3P,[0 a],[3 20]);
% % figure;imshow3(std_TE3A,[0 a],[3 20]);
% % figure;imshow3(mp_short./std_TE1A,[0 1.2],[3 20]);
% % figure;imshow3(std_TE1P./std_TE1A,[0 1.2],[3 20]);
% %% 
% dir=[pwd '\' '__20210129_224957_161000'];
% 
% seq1='SMS_STD_PART2_TE79_0007';
% path_seq1=[dir  '\' seq1];
% stdMB_TE1_ON=reshape(readDicom_SS(path_seq1),[100 100 20 4]);
% 
% seq2='SMS_STD_PART2_TE108_0008';
% path_seq2=[dir  '\' seq2];
% stdMB_TE2_ON=reshape(readDicom_SS(path_seq2),[100 100 20 4]);
% 
% seq3='SMS_STD_PART2_TE137_0010';
% path_seq3=[dir  '\' seq3];
% stdMB_TE3_ON=reshape(readDicom_SS(path_seq3),[100 100 20 4]);
% 
% seq11='STD_PART2_TE79SB_0003';
% path_seq11=[dir  '\' seq11];
% stdSB_TE1_ON=reshape(readDicom_SS(path_seq11),[100 100 20 4]);
% 
% seq22='STD_PART2_TE108SB_0004';
% path_seq22=[dir  '\' seq22];
% stdSB_TE2_ON=reshape(readDicom_SS(path_seq22),[100 100 20 4]);
% 
% seq33='STD_PART2_TE137SB_0005';
% path_seq33=[dir  '\' seq33];
% stdSB_TE3_ON=reshape(readDicom_SS(path_seq33),[100 100 20 4]);
% 
% save([directory '0129_2020_corrected']);
% %% 
% directory = 'D:\paper2\0129\'; 
% load([directory '0129_2020_corrected.mat']);
% 
% %---------TE1-------------%
% ind=4;%b=0
% figure;
% subplot(411);imshow3(stdSB_TE3_ON(:,:,:,ind),[ ],[2 10]);
% title('SB online ', 'FontSize', 10);
% subplot(412);imshow3(stdSB_P2_TE3(:,:,:,ind),[ ],[2 10]);
% title('SB offline', 'FontSize',10);
% subplot(413);imshow3(stdMB_TE3_ON(:,:,:,ind),[ ],[2 10]);
% title('MB online ', 'FontSize',  10);
% subplot(414);imshow3(stdMB_P2_TE3(:,:,:,ind),[ ],[2 10]);
% title('MB offline', 'FontSize', 10);
% 
