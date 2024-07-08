%% PRIME

% v0:       Qiang Liu, Xiaoqing Wang (Feb/05/2023), Jaejin Cho (Apr/10/2023)
% v1:       Qiang Liu (Jul/02/2023)
% v2-v9:    Yohan Jun (Aug/17/2023)

% v4 notes
%   1)  matlab/+mr/@sequence/Sequence.m
%       Line 606: inserted for separate adc events
%   2)  add r_os/r_os_low: to match b/w adcSampled/adcSampled_low

% v5 notes
%   1)  add another 3 nav acqs for 2nd echo
%   2)  modify delay TR considering 2nd echo
%   3)  add TE_low for 2nd echo
% v5_2 notes
%   1)  move 3 nav for 2nd echo to after 2nd 180 pulse

% v6 notes
%   1)  reverse 2nd echo acq for efficiency

% v7 notes
%   1)  diff 32/64 directions

% v8 notes
%   1)  frequency-encoding pF

% v9 notes
%   1)  frequency-encoding pF for 2nd echo (option)
% v9_wonav notes
%   1)  remove nav (will use dummy for LPC or DPG)

clc;  close all;  clear;

seq_file        =   'ep2d_C2_2mm_TE20_R4_2shots_slc60_dir3.seq';
save_seq_file   =   1;
seq_pns         =   1;
seq_plot        =   1;
seq_report      =   0;
check_freq      =   0;
test_check      =   1; % for debugging

bay             =   8; % Bay 2/4


%% Set system limits

B0field     =   2.89; % 6.98
lims        =   mr.opts('MaxGrad',480,'GradUnit','mT/m',...  % MaxGrad 78
                        'MaxSlew',580,'SlewUnit','T/m/s',...
                        'rfRingdownTime',10e-6,'rfDeadtime',100e-6,'B0',B0field); % max slew rate: 200 
lims1       =   mr.opts('MaxGrad',68,'GradUnit','mT/m',...  % MaxGrad 68
                        'MaxSlew',100,'SlewUnit','T/m/s',...
                        'rfRingdownTime',10e-6,'rfDeadtime',100e-6,'B0',B0field); % this lims1 is set for gz fat
lims2       =   mr.opts('MaxGrad',480,'GradUnit','mT/m',...  % MaxGrad 68
                        'MaxSlew',270,'SlewUnit','T/m/s',...
                        'rfRingdownTime',10e-6,'rfDeadtime',100e-6,'B0',B0field); % for diffusion gradients


%% Imaging parameter

% table               =   xlsread('diffusion_table/Book3_30A.xlsx'); % v2
table               =   [0 0 1; 0 1 0; 1 0 0;]; % test one axis
table               =   [0 0 0; 0 0 0; table]; % dummy + b0 + dwi % v2

% load('buda_sat_shell6_set2.mat');
% table               =   [0 0 0; table2array(buda_sat_shell6_set2)]; % dummy + dwi (b0 included) % v7
diffusion_count     =   size(table,1);

fov                 =   220e-3;
Nx                  =   110;
Nx_org              =   110;            % for gx pre-phasing calculation
Ny                  =   Nx;             % Define FOV and resolution (Nx)
thickness           =   2*1e-3;           % slice thinckness (4e-3)
Nslices             =   60;             % (28)

bval                =   1000;
bFactor             =   bval.*ones(1,size(table,1)); % s/mm^2
b0_idx              =   find(sum(table,2)==0); % v7
bFactor(b0_idx)     =   0;              % v7

% In siemens's sequence, only b0 sequence uses crushers, so in our
% sequence, we introduce a factor to keep the echo time the same for b0/b1000
crusher_switch      =   0;
crusher_d           =   0.95e-3;
% dur_rephase         =   7.6e-04;
recon               =   false;          % plot the traj for check
TE                  =   20e-3;          % (v5: 65e-3, v4: 77e-3, v2; 57e-3 ??, v1: 65e-3)
TR                  =   3800e-3;        % (v2: 3000e-3, v1: 4000e-3)

RSegment            =   2;              % multishot factor
R                   =   2.0;              % In-plane accerelation factor (2)

Echotimeshift       =   1;              % for multishot
buda                =   1;

ro_os               =   1.00;           % oversampling factor (in contrast to the product sequence we don't really need it) 1.30 (if pe_pf: 1.70) (if pe_pf/TE60: 1.10)
readoutTime         =   2.3e-4;         % this controls the readout bandwidth (v5: 8.2e-4 >> ro_os: 1.50, ro_os_low: 3.50, 8.0e-4 >> ro_os: 1.3, ro_os_low:2.6 v4: 1.1e-3, v3: 8.3e-4)

partFourierFactor           =  6/8;     % partial Fourier factor: 1: full sampling 0: start with ky=0
partFourierFactor_fe        =  6/8;     % frequency-encoding partial Fourier factor
Nx                  =   Nx * partFourierFactor_fe;

tRFex               =   3e-3;           % sec
tRFref              =   3e-3;           % sec
sat_ppm             =   -3.45;
dur_rewinder        =   1.0e-03;

deltak              =   1/fov;
deltaky             =   RSegment*R*deltak;      % Rsegement*R
kWidth              =   Nx*deltak;
kWidth_org          =   Nx_org*deltak;


%% Define Seq

seq             =   mr.Sequence(lims);      % Create a new sequence object
% rep 1 works as the dummy scan, and I also turn off the PE gradients to collect the reference
for rep =   5:diffusion_count  % 1,2 (x,y)
    if rep == 1
        pe_enable=0;           % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
    else
        pe_enable=1;
    end
    if test_check == 1
        Nmulti_end = 1;
    else
        Nmulti_end = RSegment;
    end
    for Nmulti = 1:Nmulti_end
        % fat sat
        sat_freq            =   sat_ppm*1e-6*lims.B0*lims.gamma;
        rf_fs               =   mr.makeGaussPulse(110*pi/180,'system',lims1,'Duration',8e-3,...
                                         'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
        rf_fs.phaseOffset   =   -2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
        gz_fs               =   mr.makeTrapezoid('z',lims1,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
        
        % Create 90 degree slice selection pulse and gradient
        [rf, gz, gzReph]    =   mr.makeSincPulse_QL_1(pi/2,'system',lims,'Duration',tRFex,...
                                                     'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);
        
        % Create 180 degree slice refocusing pulse and gradients
        [rf180, gz180]      =   mr.makeSincPulse(pi,'system',lims,'Duration',tRFref,...
                                                'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');
        
%         dphase              =   mr.makeDelay(dur_rephase); % we introduce this delay between the rewinder and the first diffusion gradient. % dt1_v1_public (not mentioned)
%         d_epi               =   mr.makeDelay(dur_rewinder); % 1119 % dt1_v1_public (not mentioned)

        spoiler_amp=3*8*42.58*10e2;
        est_rise=500e-6; % ramp time 280 us
        est_flat=2500e-6; %duration 600 us

        gp_r=mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gp_p=mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gp_s=mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

        gn_r=mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gn_p=mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gn_s=mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        
        % first gz180 crusher
        gz180_crusher_1     =   mr.makeTrapezoid('z',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % used for Siemens
        gz180_crusher_2     =   mr.makeTrapezoid('y',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % used for Siemens
        gz180_crusher_3     =   mr.makeTrapezoid('x',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % used for Siemens
        
        % second gz180 crusher
        gz180_c1            =   mr.makeTrapezoid('z',lims,'Amplitude',0*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref
        gz180_c2            =   mr.makeTrapezoid('y',lims,'Amplitude',-1*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref
        gz180_c3            =   mr.makeTrapezoid('x',lims,'Amplitude', 1*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref
        
        % define the output trigger to play out with every slice excitatuion
        trig                =   mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'
        
        % Phase blip in shortest possible time
        blip_dur            =   ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % round-up the duration to 2x the gradient raster time
        
        % the split code below fails if this really makes a trpezoid instead of a triangle...
        gy                  =   mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % use negative blips to save one k-space line on our way towards the k-space center
        
        % readout gradient is a truncated trapezoid with dead times at the beginnig
        % and at the end each equal to a half of blip_dur
        % the area between the blips should be defined by kWidth
        % we do a two-step calculation: we first increase the area assuming maximum
        % slewrate and then scale down the amlitude to fix the area
        extra_area          =   blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
        gx                  =   mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
        actual_area         =   gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
        gx.amplitude        =   gx.amplitude/actual_area*kWidth;
        gx.area             =   gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
        gx.flatArea         =   gx.amplitude*gx.flatTime;

        % v7 %
        extra_area_org      =   blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
%         gx_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_org+extra_area_org,'duration',readoutTime+blip_dur);
        gx_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_org+extra_area_org);
        actual_area_org     =   gx_org.area-gx_org.amplitude/gx_org.riseTime*blip_dur/2*blip_dur/2/2-gx_org.amplitude/gx_org.fallTime*blip_dur/2*blip_dur/2/2;
        gx_org.amplitude    =   gx_org.amplitude/actual_area_org*kWidth_org;
        gx_org.area         =   gx_org.amplitude*(gx_org.flatTime + gx_org.riseTime/2 + gx_org.fallTime/2);
        gx_org.flatArea     =   gx_org.amplitude*gx_org.flatTime;
        % v7 %
        
        % calculate ADC
        % we use ramp sampling, so we have to calculate the dwell time and the
        % number of samples, which are will be qite different from Nx and
        % readoutTime/Nx, respectively.
        adcDwellNyquist     =   deltak/gx.amplitude/ro_os;
        
        % round-down dwell time to 100 ns
        adcDwell            =   floor(adcDwellNyquist*1e7)*1e-7;
        adcSamples          =   floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
        
        % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
        adc                 =   mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
        
        % realign the ADC with respect to the gradient
        time_to_center      =   adc.dwell*((adcSamples-1)/2+0.5); % Siemens samples in the center of the dwell period
        adc.delay           =   round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
        % this rounding actually makes the sampling points on odd and even readouts
        % to appear misalligned. However, on the real hardware this misalignment is
        % much stronger anyways due to the grdient delays
        % adc.id              =   1; % v4
        
        % FOV positioning requires alignment to grad. raster... -> TODO
        
        % split the blip into two halves and produce a combined synthetic gradient
        gy_parts                    =   mr.splitGradientAt(gy, blip_dur/2, lims);
        [gy_blipup, gy_blipdown,~]  =   mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
        gy_blipdownup               =   mr.addGradients({gy_blipdown, gy_blipup}, lims);
        
        % pe_enable support
        gy_blipup.waveform      =   gy_blipup.waveform*pe_enable;
        gy_blipdown.waveform    =   gy_blipdown.waveform*pe_enable;
        gy_blipdownup.waveform  =   gy_blipdownup.waveform*pe_enable;
        
        % phase encoding and partial Fourier
        Ny_pre      =   round((partFourierFactor-1/2)*Ny-1);  % PE steps prior to ky=0, excluding the central line
        Ny_pre      =   round(Ny_pre/RSegment/R);
        Ny_post     =   round(Ny/2+1); % PE lines after the k-space center including the central line
        Ny_post     =   round(Ny_post/RSegment/R);
        Ny_meas     =   Ny_pre+Ny_post;
        
        % pre-phasing gradients
        gxPre           = mr.makeTrapezoid('x',lims,'Area',-gx_org.area/2); % v7
        gyPre           = mr.makeTrapezoid('y',lims,'Area',(Ny_pre*deltaky-(Nmulti-1)*R*deltak));
        [gxPre,gyPre]   = mr.align('right',gxPre,'left',gyPre);
        % relax the PE prepahser to reduce stimulation
        gyPre           = mr.makeTrapezoid('y',lims,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
        % gxPre's duration is always large because we use PF for PE, so for
        % ms EPI we are safe -- we don't need to notice this duration changing --QL
        gyPre.amplitude = gyPre.amplitude*pe_enable;
        
        % post gradient to go back to k-space center
        gyPre_post_area = -1*(Ny_post-0.5)*deltaky; % this seems to be correct
        
        if(mod(Nmulti,2) ~= 0)
            gyPre_post_area = (Ny_post-1)*deltaky; % this seems to be correct
        end
        
        gxPre_post_area_nav = -1*(0.5*gx.flatArea + gx.amplitude*gx.fallTime/2); % partFourierFactor should be larger than 0.5 % v7
        gxPre_post_area     = -1*(abs(gxPre_post_area_nav*2) - (0.5*gx_org.flatArea + gx_org.amplitude*gx_org.fallTime/2)); % partFourierFactor should be larger than 0.5 % v7
        gxPre_post_area_nav = -1*(abs(gxPre_post_area_nav*2) - (0.5*gx_org.flatArea + gx_org.amplitude*gx_org.fallTime/2)); % partFourierFactor should be larger than 0.5 % v7
        
        if(mod(Ny_meas,2) == 0)
%             gxPre_post_area_nav = 0.5*gx.flatArea + gx.amplitude*gx.riseTime/2; % v7
            gxPre_post_area     = 0.5*gx_org.flatArea + gx_org.amplitude*gx_org.riseTime/2; % v7
        end
        
%         gxPre_post_nav  = mr.makeTrapezoid('x',lims,'Area', gxPre_post_area_nav, 'Duration', 0.0012); % Fixme: gxPre_post_area is still not centered % v7 % v9_wonav_uc
        gxPre_post      = mr.makeTrapezoid('x',lims,'Area', gxPre_post_area, 'Duration', 0.0012); % Fixme: gxPre_post_area is still not centered % v7
        gyPre_post      = mr.makeTrapezoid('y',lims,'Area', gyPre_post_area, 'Duration', 0.0012); % (QL, XW)
        
        [gxPre_post,gyPre_post]=mr.align('right',gxPre_post,'left',gyPre_post);
        % relax the PE prepahser to reduce stimulation
        gyPre_post = mr.makeTrapezoid('y',lims,'Area',gyPre_post.area,'Duration',mr.calcDuration(gxPre_post,gyPre_post));
        
        gyPre_post.amplitude=gyPre_post.amplitude*pe_enable;
        
        % BUDA
        if (buda)
            if (mod(Nmulti,2) == 0)
                gyPre           =   mr.scaleGrad(gyPre,-1);
                gy_blipdownup   =   mr.scaleGrad(gy_blipdownup,-1);
                gy_blipdown     =   mr.scaleGrad(gy_blipdown,-1);
                gy_blipup       =   mr.scaleGrad(gy_blipup,-1);
            end
        end
        
        dur2end=(Ny_post-0.5)*mr.calcDuration(gx)+mr.calcDuration(gxPre_post,gyPre_post);

        %% Calculate delay times, QL
        
        if (rep < 3) % need to adjust delta TE1 and delta TE2 for b0, bcuz of the crusher
            crusher_switch=1;
        else
            crusher_switch=0;
        end
        
        % Calculate delay times
        durationToCenter        = (Ny_pre + 0.5) * mr.calcDuration(gx);
        durationToCenter2       = (Ny_post + 0.5) * mr.calcDuration(gx); % v4
        
        rfCenterInclDelay       = rf.delay + mr.calcRfCenter(rf);
        rf180centerInclDelay    = rf180.delay + mr.calcRfCenter(rf180);
        
        delayTE1                = ceil((TE/2 - mr.calcDuration(rf,gz) - mr.calcDuration(gzReph) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
%         delayTE1                = delayTE1 - mr.calcDuration(gxPre) - mr.calcDuration(gxPre_post_nav) - 3*mr.calcDuration(gx); % JC added for Nav % v9_wonav_uc
%         delayTE1                = delayTE1 - mr.calcDuration(gxPre_low) - mr.calcDuration(gxPre_low_post) - 3*mr.calcDuration(gx_low); % v5 Nav_low
        delayTE1                = delayTE1 - crusher_switch * mr.calcDuration(gz180_crusher_1); %QL
%         delayTE1                = delayTE1 - dur_rephase; % QL (not mentioned in dt1_v1_public, but needed for accurate TE/2)
        assert(delayTE1>=0);
        
        gxPre.delay=0;
        gyPre.delay=0;
        delayTE2                = ceil((TE/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
        delayTE2                = delayTE2 - mr.calcDuration(gxPre,gyPre);
        delayTE2                = delayTE2 - crusher_switch*mr.calcDuration(gz180_crusher_1); %QL
        assert(delayTE2>=0);
                
        [gxPre,gyPre]           = mr.align('right',gxPre,'left',gyPre);
        
        % diffusion weithting calculation
        % delayTE2 is our window for small_delta
        % delayTE1+delayTE2-delayTE2 is our big delta
        % we anticipate that we will use the maximum gradient amplitude, so we need
        % to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
        small_delta=delayTE2-ceil(lims2.maxGrad/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        big_delta=delayTE1+mr.calcDuration(rf180,gz180);
        % we define bFactCalc function below to eventually calculate time-optimal
        % gradients. for now we just abuse it with g=1 to give us the coefficient
        % g=sqrt(bFactor*1e6/bFactCalc(1,small_delta,big_delta)); % for now it looks too large!
        g=sqrt(bFactor(1,rep)*1e6/bFactCalc(1,small_delta,big_delta)); % QL

        gr=ceil(g/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims2); % v2
        

        %% split the diffusion gradient to 3 axis, new version based on Jon and Maxim's suggestions (QL) % v2

        g_x=g.*table(rep,1);
        g_y=g.*table(rep,2);
        g_z=g.*table(rep,3);

        if ((sum(table(rep,:))==0)||(sum(table(rep,:))==1)) % b=0 or dwi with diffusion gradient on one axis we keep using the older version
            g_xr=ceil(abs(g_x)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime; % QL: diffusion Oct 27
            g_yr=ceil(abs(g_y)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            g_zr=ceil(abs(g_z)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            gDiff_x=mr.makeTrapezoid('x','amplitude',g_x,'riseTime',g_xr,'flatTime',small_delta-g_xr,'system',lims2);
            gDiff_y=mr.makeTrapezoid('y','amplitude',g_y,'riseTime',g_yr,'flatTime',small_delta-g_yr,'system',lims2);
            gDiff_z=mr.makeTrapezoid('z','amplitude',g_z,'riseTime',g_zr,'flatTime',small_delta-g_zr,'system',lims2);
            duration_diff=max(max(mr.calcDuration(gDiff_x), mr.calcDuration(gDiff_y)), mr.calcDuration(gDiff_z));

        else % 3 axis, 'rotation'
            [azimuth,elevation,r] = cart2sph(g_x,g_y,g_z);
            polar= -(pi/2-elevation);

            Gr=mr.rotate('z',azimuth,mr.rotate('y',polar,gDiff));
            if size(Gr,2)==3
                gDiff_x=Gr{1,2};
                gDiff_y=Gr{1,3};
                gDiff_z=Gr{1,1};
            else
                if size(Gr,2)==2
                    diffusion_blank=find( table(rep,:)==0);
                    switch diffusion_blank
                        case 2
                            gDiff_x=Gr{1,2};
                            gDiff_z=Gr{1,1};
                            gDiff_y=gDiff; gDiff_y.channel='y'; gDiff_y.amplitude=0; gDiff_y.area=0; gDiff_y.flatArea=0;
                        case 1
                            gDiff_z=Gr{1,1};
                            gDiff_y=Gr{1,2};
                            gDiff_x=gDiff; gDiff_x.amplitude=0; gDiff_x.area=0; gDiff_x.flatArea=0;gDiff_x.channel='x';
                        case 3
                            gDiff_x=Gr{1,2};
                            gDiff_y=Gr{1,1};
                            gDiff_z=gDiff; gDiff_z.amplitude=0; gDiff_z.area=0; gDiff_z.flatArea=0;gDiff_z.channel='z';
                    end
                end
            end
            duration_diff=mr.calcDuration(gDiff);
        end

        assert(duration_diff<=delayTE1);
        assert(duration_diff<=delayTE2);
        
        
        %% Calculate the echo time shift for multishot EPI (QL)
        actual_esp=gx.riseTime+gx.flatTime+gx.fallTime;
        TEShift=actual_esp/RSegment;
        TEShift=round(TEShift,5); % from Berkin: roundn didn't work for the latest matlab, changed to round (sign -/+) % v2
        TEShift_before_echo=(Nmulti-1)*TEShift;
        if TEShift_before_echo ==0
            TEShift_before_echo=0.00001; % apply the minimum duration for the no delay case
        end
        TEShift_after_echo=(RSegment-(Nmulti-1))*TEShift;
        dETS_before=mr.makeDelay(TEShift_before_echo);
        dETS_after=mr.makeDelay(TEShift_after_echo);
        
        
        %% TR calculation
        
        delayTR     =   ceil((TR/Nslices - mr.calcDuration(gp_r) - mr.calcDuration(gn_r) - TE - rfCenterInclDelay-dur2end)/seq.gradRasterTime)*seq.gradRasterTime; % v5
%         delayTR     =   delayTR/Nslices;
        delayTR     =   round(delayTR,3); % pulseq will allow delay time at power -3 % v2
        dTR         =   mr.makeDelay(delayTR);
        
        
        %% interleaved

        for islice_1=1:Nslices
            freqOffset_factor(islice_1)=islice_1-1-(Nslices-1)/2;
        end
        slic_indexS=[2:2:Nslices 1:2:Nslices];
        interleaved_freqOffset_factor=freqOffset_factor(slic_indexS);
        
        
        %% Define Seq Block
        
        if test_check == 1
            slc_end = 1;
        else
            slc_end = Nslices;
        end
        for slc = 1:slc_end
            
            % 90 RF Module Define
            seq.addBlock(gp_r,gp_p,gp_s); % v11
            seq.addBlock(rf_fs,gn_r,gn_p,gn_s); % v11

            rf.freqOffset           =   gz.amplitude*thickness*interleaved_freqOffset_factor(slc);
            %     rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
            rf180.freqOffset        =   gz180.amplitude*thickness*interleaved_freqOffset_factor(slc);
            %     rf180.phaseOffset=pi/2-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
            rf180.phaseOffset       =   pi/2;
            seq.addBlock(rf,gz,trig);
            seq.addBlock(gzReph); %rewinder QL
%             seq.addBlock(dphase); % this delay split the rewinder and the diffusion gradient (v2: not mentioned in v1)
            
            % nav (JC)
%             seq.addBlock(gxPre); % time need to be recalculated % v9_wonav_uc
%             seq.addBlock(gx,adc); % v9_wonav_uc
%             gx.amplitude            =   -gx.amplitude; % v9_wonav_uc
%             seq.addBlock(gx,adc); % v9_wonav_uc
%             gx.amplitude            =   -gx.amplitude; % v9_wonav_uc
%             seq.addBlock(gx,adc); % v9_wonav_uc
%             seq.addBlock(gxPre_post_nav); % v9_wonav_uc
            
            % diff gradient
            seq.addBlock(mr.makeDelay(delayTE1),gDiff_x,gDiff_y,gDiff_z);
            
            % 180 RF Module Define
            if (rep < 3)
                if(~recon)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
            end
            
            seq.addBlock(rf180,gz180); % QL: it used to be gz180n
            
            if (rep < 3)
                if(~recon)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
            end
            
            % diff gradient
            seq.addBlock(mr.makeDelay(delayTE2),gDiff_x, gDiff_y, gDiff_z);
            
            % QL for no blip
            if rep == 1
                gy_blipup.last=0;
                gy_blipdown.first=0;
                gy_blipdownup.last=0;
                gy_blipdownup.first=0;
            end
            
            if (Echotimeshift)
                seq.addBlock(dETS_before); % echotimeshift for multishot
            end
            
            seq.addBlock(gxPre,gyPre);
            
            for i=1:Ny_meas
                if i==1
                    seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_meas
                    seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                else
                    seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
                % gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                % QL: for odd pe lines, the k-space trajectory was not correct
                if (mod(Ny_meas,2)==1)
                    if (i<Ny_meas)
                        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    end
                else
                    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                end
            end
            
            seq.addBlock(gxPre_post,gyPre_post); % go back to k-space center
            
            if (Echotimeshift)
                seq.addBlock(dETS_after); % echotimeshift for multishot
            end
            
            seq.addBlock(dTR); % seperate
        end %slice loop
    end % for multishot
    disp(['diffusion dir:', num2str(rep)])
    clear gDiff gDiff_x gDiff_y gDiff_z
end % EPI REF diffusion LOOP


%% check whether the timing of the sequence is correct

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end


%% do some visualizations

if seq_plot == 1
    tic
    fprintf('Plotting figures... ');

    seq.plot(); % plot sequence waveforms
    
    % trajectory calculation
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
    
    % plot k-spaces
    figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
    hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
    
    figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    hold on; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
    toc
end


%% prepare the sequence output for the scanner

seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi-diff');

if seq_pns == 1
    tic
    fprintf('Checking PNS... ');
    % pns < 80\% (human)
    [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/MP_GradSys_P034_c_CX600.asc'); % PRISMA-XR
    pns = max(pns_n(:));
    
    clear ktraj_adc t_adc ktraj t_ktraj t_excitation t_refocusing slicepos t_slicepos pns_n pns_c tpns
    toc
end

% if pns_ok
    if save_seq_file == 1
        tic
        fprintf('Saving .seq file... ');
        seq.write(seq_file);
        save(strcat(seq_file(1:end-4),'_param'));
        toc
    end
% end

% seq.install('siemens');
% seq.sound(); % simulate the seq's tone


%% check forbidden frequencies

if check_freq == 1
    tic
    fprintf('Checking frequencies... ');
    sys = lims;
    gradSpectrum;
    toc
end


%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

if seq_report == 1
    tic
    fprintf('Checking test report... ');
    rep = seq.testReport;
    fprintf([rep{:}]);
    toc
end


%%

function b=bFactCalc(g, delta, DELTA)
% see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 39????5 (2012) DOI 10.1002/cmr.a
% b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
% in pulseq we don't need gamma as our gradinets are Hz/m
% however, we do need 2pi as diffusion equations are all based on phase
% for rect gradients: sigma=1 lambda=1/2 kappa=1/3
% for trapezoid gradients: TODO
sigma=1;
%lambda=1/2;
%kappa=1/3;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end