% add Pulseq, Toppe, and PulCeq toolbox
% for Toppe interpreter TV6
% From Jon

system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
system('git clone git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
system('git clone git@github.com:toppeMRI/toppe.git');
addpath toppe