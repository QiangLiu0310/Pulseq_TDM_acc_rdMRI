function simse
% Create 90-180 spin-echo pair and simulate

slThick = 0.24;  % cm
sysGE = toppe.systemspecs();

% 90 excitation pulse
tbw = 4;
dur = 4;        % ms
[rfex, gzex] = toppe.utils.rf.makeslr(90, slThick, tbw, dur, eps, sysGE, 'type', 'ex', 'writeModFile', false);

% refocusing pulse and crusher
[rfr, gzr] = toppe.utils.rf.makeslr(180, 1.2*slThick, tbw, dur, eps, sysGE, 'type', 'se', 'writeModFile', false);
ncycles = 2; 
gcrush = toppe.utils.makecrusher(ncycles,slThick,sysGE); 

% apply slice offset
sliceOffset = -7.2;  % cm
dt = 4e-6;  % s
f = sysGE.gamma*max(gzex)*1e-4*sliceOffset;   % Hz
t = dt*(1:length(rfex))';
rfex = rfex.*exp(1i*2*pi*f*t);
iRfCenter = find(abs(rfex)==max(abs(rfex)))-1;
rfex = rfex*exp(-1i*angle(rfex(iRfCenter)));

f = sysGE.gamma*max(gzr)*1e-4*sliceOffset;   % Hz
t = dt*(1:length(rfr))';
rfr = rfr.*exp(1i*2*pi*f*t);
iRfCenter = find(abs(rfr)==max(abs(rfr)));
rfr = rfr*exp(-1i*angle(rfr(iRfCenter)));

% apply pi/2 phase offset for excitation pulse since blocksim rotates by pi/2
rfex = rfex*exp(1i*pi/2);

% apply 'wrong' phase offset for refocusing pulse to check effect on slice profile
ph_ref = pi/4;
rfr = rfr*exp(1i*ph_ref);

% simulation parameters
fov = 1;            % cm
m0 = [0 0 1];       % initial magnetization
z = sliceOffset + linspace(-fov/2, fov/2, 500);
T1 = inf;           % ms
T2 = inf;            % ms
dt = 4e-3;           % ms

% simulate 90 by itself first (to check phase)
figure;
rf = [rfex];
gz = [gzex];
[m] = toppe.utils.rf.slicesim(m0, rf, gz, dt, z, T1, T2, true);
ttl = sprintf('90-zoff-%.1fcm', sliceOffset);
subplot(131); title(ttl);
print([ttl '.png'], '-dpng');

% simulate 90-180 pair with crushers
figure;
rf = [rfex; 0*gcrush; rfr; 0*gcrush];
gz = [gzex; gcrush; gzr; gcrush];
[m] = toppe.utils.rf.slicesim(m0, rf, gz, dt, z, T1, T2, true);
ttl = sprintf('90-180-zOffset%.1fcm-phsOffset%.2f', sliceOffset, ph_ref);
subplot(131); title(ttl);
print([ttl '.png'], '-dpng');

% simulate 90-180 pair with b1 ~= 1.0
figure;
b1_inhomo = 0.8;
[m] = toppe.utils.rf.slicesim(m0, b1_inhomo*rf, gz, dt, z, T1, T2, true);
ttl = sprintf('90-180-zOffset%.1fcm-phsOffset%.2f-b1inhomo-%.2f', sliceOffset, ph_ref, b1_inhomo);
subplot(131); title(ttl);
print([ttl '.png'], '-dpng');

