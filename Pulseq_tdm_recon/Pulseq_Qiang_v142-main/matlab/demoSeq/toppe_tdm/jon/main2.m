%% Get MATLAB toolboxes

% Pulseq (+mr namespace)
%system('git clone --branch v1.4.2 git@github.com:pulseq/pulseq.git');
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% tools for converting Pulseq file to GE scan files
system('git clone --branch main git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
system('git clone --branch main git@github.com:toppeMRI/toppe.git');
addpath toppe


%% create .seq file
%makeSpiralSequence2;
%writeEpiDiffusionRS_3shot_v3_epi_in_1_block_appa_magnus;
writeEpiDiffusionRS_3shot_v2_2mm_appa_v2_jfn;


%% convert to .tar file for GE
fn = 'tdm_3e_mb1_21sli_2mm_R3_test';
fn = 'epidiff_3_shot_ref_2mm_21sli_1';
% ceq = seq2ceq('mp2ragesagte.seq', 'ignoreSegmentLabels', true);
ceq = seq2ceq([fn '.seq']);
% save ceq ceq
% load ceq

sysGE = toppe.systemspecs('rfDeadTime', 0, ...
                          'rfRingdownTime', 0, ...
                          'adcDeadTime', 0, ...
                          'maxGrad', 30, 'maxSlew', '80', ...
                          'maxRF', 0.25, ...
                          'psd_rf_wait', 0, ...
                          'psd_grd_wait', 0);

%ceq2ge(ceq, sysGE, [fn '.tar']);
ceq2ge(ceq, sysGE, 'out.tar');

system('tar xf out.tar');

