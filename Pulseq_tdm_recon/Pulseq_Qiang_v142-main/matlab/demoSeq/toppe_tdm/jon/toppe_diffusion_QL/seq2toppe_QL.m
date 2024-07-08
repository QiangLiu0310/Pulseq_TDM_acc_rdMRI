clear;close all;clc;
%% convert to .tar file for GE
fn = 'signleecho_mb1_21sli_2mm_R3_te3';
ceq = seq2ceq([fn '.seq']);

sysGE = toppe.systemspecs('rfDeadTime', 0, ...
                          'rfRingdownTime', 60, ... % QL: 0
                          'adcDeadTime', 0, ...
                          'maxGrad', 30, 'maxSlew', '80', ...
                          'maxRF', 0.25, ...
                          'psd_rf_wait', 60, ...
                          'psd_grd_wait', 0);


ceq2ge(ceq, sysGE, 'scan6_2.tar');

% system('tar xf scan4.tar');

