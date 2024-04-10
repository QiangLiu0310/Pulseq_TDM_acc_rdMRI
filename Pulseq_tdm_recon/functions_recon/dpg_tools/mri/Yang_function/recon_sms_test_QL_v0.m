% test for XA30
% Qiang Liu based on Berkin's code

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/data_processing/read_meas_dat__20140924112147'))
addpath(genpath('/data/pnl/home/ql087/data_processing/FID-A-master'))


file_path = '/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_15_bwh_prisma_phantom/raw/';
file_name='meas_MID00863_FID23317_Diffusion_SMS2_R2_PA.dat';


tic
    dat = mapVBVD([file_path, file_name]);
toc

tic
    meas.prot = read_meas_prot_bb_v0([file_path, file_name]);
toc

prot = meas.prot;

if isempty(prot.aflRegridADCDuration)
    prot.aflRegridADCDuration = dat{2}.hdr.Meas.aflRegridADCDuration(1);
end

[ image ]=recon_sms_std_xa30_0(dat,prot);
