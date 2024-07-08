%% Convert to TOPPE .tar file for execution on GE scanners
fn = 'mprage.seq';

% First convert to Matlab struct based on PulCeq
ceq = seq2ceq(fn);

% Get sys struct from .seq file (Ideally remove this)
seq = mr.Sequence(); seq.read(fn); sys = seq.sys;

% Geometry
% FOV = [19.2 25.6 25.6].*1e-2; % Field of view (m)
% res = [1 1 1].*1e-3;          % Resolution (m)
% N = FOV./res;                 % Matrix size

% TOPPE system specs struct
maxGrad = sys.maxGrad/sys.gamma*100;   % G/cm
maxSlew = 1.01*sys.maxSlew/sys.gamma/10;    % G/cm/ms. Factor > 1 is fudge factor to avoid exceeding limit after interpolating to GE raster time.
sysGE = toppe.systemspecs('maxGrad', maxGrad, ...   % G/cm
    'maxSlew', maxSlew, ...                         % G/cm/ms
    'maxRF', 0.15, ...
    'maxView', N(1), ...                            % Length of second fastest dimension
    'rfDeadTime', sys.rfDeadTime*1e6, ...           % us
    'rfRingdownTime', sys.rfRingdownTime*1e6, ...   % us
    'adcDeadTime', sys.adcDeadTime*1e6);            % us

% Convert to GE
ceq2ge(ceq,sysGE,'mprage.tar','verbose',true);
