%% Creates a 3D MP-RAGE scan using pulseq with optional acceleration via GRAPPA
% Original code by Berkin Bilgic
% Generalized/repurposed by Rex Fung

%% User parameters (change these for your need

% Select scanner ('Siemens' or 'GE')
scanner = 'GE';

% Geometry
FOV = [19.2 25.6 25.6].*1e-2; % Field of view (m)
res = [1 1 1].*1e-3;          % Resolution (m)
N = FOV./res;                 % Matrix size

% Acceleration
R = 1; % Undersampling factor in slowest (slice) dimension

% Test report
% optional slow step, but useful for testing during development
% e.g. for the real TE, TR or for staying within slew rate limits  
run_test = false;

%% Experimental parameters
% Scanner hardware limits
% Siemens 3T: rfDeadTime = 100e-6; rfRingdownTime = 10e-6; adcDeadTime = 10e-6;
% GE 3T: rfDeadTime = 50e-6; rfRingdownTime = 60e-6; adcDeadTime = 40e-6;
if strcmp(scanner,'Siemens')
    sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
        'MaxSlew', 70, 'SlewUnit', 'T/m/s', ...
        'B0', 123/128*3, ...
        'rfDeadTime', 100e-6, ...
        'rfRingdownTime', 10e-6, ...
        'adcDeadTime', 10e-6);
elseif strcmp(scanner,'GE')
    sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
        'MaxSlew', 70, 'SlewUnit', 'T/m/s', ...
        'B0', 3.0, ...
        'rfDeadTime', 50e-6, ...
        'rfRingdownTime', 60e-6, ...
        'adcDeadTime', 40e-6);
else
    fprintf('Please select scanner')
    return;
end
seq = mr.Sequence(sys);   % Create a new sequence object

% Common params
alpha = 8;                % Flip angle
dwell = 16e-6;            % ADC sample time (s). For GE, must be multiple of 2us

% Additional params
ro_spoil = 3;             % Additional k-max excursion for RF spoiling
TI = 1.07;                % Inversion time (s). To match the ABCD protocol. Was: 1.1
TRout = 2.5;              % Outer loop repetition time (s)
rfSpoilingInc = 117;      % RF spoiling phase increment (degrees)

% Axes params
ax=struct;                % encoding axes
ax.d1='z';                % the fastest dimension (readout)
ax.d2='x';                % the second-fastest dimension (the inner pe loop)
ax.d3=setdiff('xyz',[ax.d1 ax.d2]); % automatically set the slowest dimension
ax.n1=strfind('xyz',ax.d1);
ax.n2=strfind('xyz',ax.d2);
ax.n3=strfind('xyz',ax.d3);

% Readout params
ro_dur = N(ax.n1) * dwell; % readout duration
ro_os = 1;                 % readout oversampling rate

%% Create alpha-degree hard or binominal water excitation pulse
% Excitation params
fatChemShift = 3.5; % ppm
gamma = 4.2576e7;   % Hz/Tesla
gamG = gamma/1e4;    % Hz/Gauss
fatOffres = gamma*sys.B0*fatChemShift*1e-6; % fat resonance frequency offset (Hz)

% Excitation struct
ex.mode = 'hard';  % 'we' = 1-1 binomial water excitation; 'hard' = rectangular pulse
ex.nhard = 50; % number of waveform samples in each alpha/2 hard pulse
ex.hard = (alpha/2/360) / (gamG * ex.nhard * sys.rfRasterTime) * ones(ex.nhard,1);
if strcmp(ex.mode, 'we')
    ex.ngap = round(1/fatOffres/2/sys.rfRasterTime - ex.nhard); 
else
    ex.ngap = 0;
end
ex.signal = [ex.hard; zeros(ex.ngap,1); ex.hard];

% Magnetization pulse params
tmp = length(ex.signal) * sys.rfRasterTime;  % sec
tpad = sys.gradRasterTime - mod(tmp, sys.gradRasterTime);
npad = round(tpad/sys.rfRasterTime);
ex.signal = [ex.signal; zeros(npad,1)]; % pad to Siemens gradient/block duration raster time
n = length(ex.signal);  % temporary shorthand
rf = mr.makeArbitraryRf(ex.signal, alpha/180*pi, 'system', sys);

%% Create inversion pulse
% TODO: make adiabatic pulse
% rf180 = mr.makeAdiabaticPulse('hypsec', sys, 'Duration', 10.24e-3, 'dwell', 1e-5);
n180 = 1000;
rf180 = 180/360 / (gamG * n180 * sys.rfRasterTime) * ones(n180, 1);
rf180 = mr.makeArbitraryRf(rf180, pi, 'system', sys);

%% Define other gradients and ADC events
deltak=1./FOV;
gro = mr.makeTrapezoid(ax.d1, ...
    'Amplitude', N(ax.n1)*deltak(ax.n1)/ro_dur, ...
    'FlatTime', ceil(ro_dur/sys.gradRasterTime)*sys.gradRasterTime, ...
    'system',sys);
adc = mr.makeAdc(N(ax.n1)*ro_os, ...
    'Duration', ro_dur,...
    'Delay', gro.riseTime, ...  
    'system',sys);
groPre = mr.makeTrapezoid(ax.d1, ...
    'Area', -gro.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gro.riseTime),...
    'system',sys); % the first 0.5 is necessary to acount for the Siemens sampling in the center of the dwell periods
gpe1 = mr.makeTrapezoid(ax.d2, ...   
    'Area', deltak(ax.n2)*(N(ax.n2)/2),...  
    'system',sys);  % maximum PE1 gradient
gpe2 = mr.makeTrapezoid(ax.d3, ...
    'Area', deltak(ax.n3)*(N(ax.n3)/2), ...
    'system',sys);  % maximum PE2 gradient
gslSp = mr.makeTrapezoid(ax.d3, ...
    'Area', max(deltak.*N)*4, ...  % spoil with 4x cycles per voxel
    'Duration', 10e-3, ...
    'system',sys);  

% we cut the RO gradient into two parts for the optimal spoiler timing
[gro1,groSp]=mr.splitGradientAt(gro,gro.riseTime+gro.flatTime);
% gradient spoiling
if ro_spoil>0
    groSp = mr.makeExtendedTrapezoidArea(gro.channel, gro.amplitude, 0, deltak(ax.n1)/2*N(ax.n1)*ro_spoil, sys);
end

% calculate timing of the fast loop 
% we will have two blocks in the inner loop:
% 1: RF 
% 2: prewinder,phase enconding + readout + spoilers/rewinders
[groPre,~] = mr.align('right',groPre,mr.makeDelay(mr.calcDuration(gpe1,gpe2)-gro.riseTime));
gro1.delay = mr.calcDuration(groPre);
groSp.delay = mr.calcDuration(gro1);
adc.delay = gro1.delay+gro.riseTime;
gro1 = mr.addGradients({gro1,groPre,groSp},'system',sys);
gpe1c = mr.addGradients({gpe1, ...
    mr.makeTrapezoid(ax.d2, 'Area', -gpe1.area, 'duration', groSp.shape_dur, 'delay', groSp.delay, 'system',sys)});
gpe2c = mr.addGradients({gpe2, ...
    mr.makeTrapezoid(ax.d3, 'Area', -gpe2.area, 'duration', groSp.shape_dur, 'delay', groSp.delay, 'system',sys)});
TRinner = mr.calcDuration(rf)+mr.calcDuration(gro1); % we'll need it for the TI delay

% peSteps -- control reordering
pe1Steps = ((0:N(ax.n2)-1)-N(ax.n2)/2)/N(ax.n2)*2;
pe2Steps = ((0:N(ax.n3)-1)-N(ax.n3)/2)/N(ax.n3)*2;

% TI calc
TIdelay=round((TI-(find(pe1Steps==0)-1)*TRinner-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay)-rf.delay-mr.calcRfCenter(rf))/sys.blockDurationRaster)*sys.blockDurationRaster;
TRoutDelay=TRout-TRinner*N(ax.n2)-TIdelay-mr.calcDuration(rf180);

% pre-register objects that do not change while looping
gslSp.id=seq.registerGradEvent(gslSp);
gro1.id=seq.registerGradEvent(gro1);
[~, gpe1c.shapeIDs]=seq.registerGradEvent(gpe1c);
[~, gpe2c.shapeIDs]=seq.registerGradEvent(gpe2c);
[~, rf.shapeIDs]=seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes 
[rf180.id, rf180.shapeIDs]=seq.registerRfEvent(rf180); % 

% define GRAPPA undersampling pattern along slowest dimension (outer pe loop)
nACS = N(ax.n3)/8;  % number of auto-calibration lines (dense sampling in center)
S = 0*(1:N(ax.n3));
S(1:R:end) = 1;
S( (end/2-nACS/2):(end/2+nACS/2-1) ) = 1;
J = find(S == 1);

% Loop over phase encodes and define sequence blocks
% j < 0: Dummy shots to reach steady state
% j = 0: ADC is turned on and used for receive gain calibration on GE scanners
% j > 0: Image acquisition
nDummyShots = 2;

% Scan loop
fprintf('\nWriting scanloop: pe line 0 of %d', N(ax.n3))
prev_j = 1; % tracker for clearing previous output
for j = [-nDummyShots:0 J] % Slice loop (y)

    % Progress message
    for ib = 1:strlength(sprintf('Writing scanloop: pe line %d of %d', prev_j, N(ax.n3)))
        fprintf('\b');
    end
    prev_j = j;
    fprintf('Writing scanloop: pe line %d of %d', j, N(ax.n3))

    % inversion pulse
    segmentID = 1; % Same for all j
    seq.addBlock(rf180, mr.makeLabel('SET', 'TRID', segmentID));

    % spoiler and delay
    seq.addBlock(mr.makeDelay(TIdelay), gslSp);

    % RF spoiling
    rf_phase = 0;
    rf_inc = 0;   

    % pre-register the PE gradients that repeat in the inner loop
    if j > 0
        gpe2cj=mr.scaleGrad(gpe2c,pe2Steps(j));
        gpe2cj.id=seq.registerGradEvent(gpe2cj);
    end

    segmentID = 2 + 3*((j >= 0) + (j > 0)); % Different segment ID for dummy, calibration, and acquisition shots
    for i = 1:N(ax.n2) % PE loop (x)
        % excitation
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        seq.addBlock(rf, mr.makeLabel('SET', 'TRID', segmentID));

        % readout (turn off ADC for dummy scans, and PE gradient for prescan)
        if j < 0
            seq.addBlock(gro1);
        elseif j == 0
            seq.addBlock(adc,gro1);
        else
            seq.addBlock(adc,gro1,mr.scaleGrad(gpe1c,pe1Steps(i)),gpe2cj);
        end
        
        % RF spoiling
        rf_inc=mod(rf_inc + rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase + rf_inc, 360.0);
    end
    
    % TR delay
    segmentID = 3; % Same for all j
    seq.addBlock(mr.makeDelay(TRoutDelay), mr.makeLabel('SET', 'TRID', segmentID));
end

% Add noise scans.
% Do this last, since receive gain for GE is based on signal from 
% beginning of sequence.
% First insert ~5s pause to allow magnetization/system to settle.
% Set delay of adc block so duration is multiple of sys.gradRasterTime.
if false
    adcDur = adc.numSamples * adc.dwell; % sec
    adcDurNew = adcDur + sys.gradRasterTime - mod(adcDur, sys.gradRasterTime);
    adc.delay = adc.delay + (adcDurNew - adcDur);
    seq.addBlock(mr.makeDelay(5)); % sec
    for i = 1:N(ax.n2)  
        seq.addBlock(adc);
    end
end

fprintf('. Sequence ready\n');


%% check sequence timing
fprintf('Checking Pulseq timing... ');
[ok, error_report]=seq.checkTiming();

if (ok)
    fprintf('Pulseq timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% plot, etc
seq.plot('TimeRange',[0 TRout*5]); 

%%
seq.setDefinition('FOV', FOV);
seq.setDefinition('Name', 'mprage');

seq.write('mprage.seq') % Write to pulseq file
% seq.install('siemens');

%% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
if max(N) <= 32
    tic;
    [kfa,ta,kf]=seq.calculateKspacePP();
    toc
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
end

%% Test report
if run_test
    rep = seq.testReport; 
    fprintf([rep{:}]);
end
