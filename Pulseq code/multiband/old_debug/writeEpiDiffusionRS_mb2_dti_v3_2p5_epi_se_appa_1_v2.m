% undersampling R2 version 1 Qiang Liu 08042022
% then we try the multishot version based on R2 version 08102022
% echo time shift is used in ms-EPI we just followed the way adding a delay
% 0-1/R ESP-2/R ESP... we seperate the delay into Nsegments parts, and we
% can multiply it by Nsegment when apply it
% This version contains: multishot+in-plane accerelation+dti
% Qiang Liu Feb/5/2023
% add PA for b0
% add a MB2 version Dec/25/2023
% split into 3 sub-groups a, b, and c

clc;close all;clear all;

seq_file='signleecho_mb2_60sli_2p5mm_R3_multib_te1_a_v3.seq'; % te [1, 2, 3], group a, b
% Set system limits
B0field=2.89;% 6.98
lims = mr.opts('MaxGrad',58,'GradUnit','mT/m',...
    'MaxSlew',180,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0',B0field);  %max slew rate: 200
lims1 = mr.opts('MaxGrad',58,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', B0field);  % this lims1 is set for gz fat
lims2 = mr.opts('MaxGrad',78,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', B0field);  % for diffusion gradients

tmp=load('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Pulseq_Qiang_PNL/matlab/demoSeq/diffusion_table/TMD_part1_bvec.txt');
tmp=tmp(3:end,:);
tmp(:,3)=-tmp(:,3);
table=zeros(25,3,3); % dir-xyz-bvalue
for i_bval=1:3
    table(:,:,i_bval)=tmp((i_bval:3:end),:);
end
clear i_bval tmp

table=[zeros(3,3,3); table];
crusher_switch=0;
bFactor=[1000:1000:3000]'.*ones(1,size(table,1)); % s/mm^2
bFactor(:,1:3)=0;
dummy_count=size(table,1);

seq=mr.Sequence(lims);     % Create a new sequence object
fov=220e-3; Nx=88; Ny=Nx; % Define FOV and resolution
thickness=2.5e-3;            % slice thinckness
Nslices=60;
TE1=70e-3;
TE2=95e-3;
TE3=120e-3;
TE=TE1;
TE_index=1; % 1,2,3
mb=2;
cai_flag=1;
RSegment=1;% multishot factor
R=3;% In-plane accerelation factor
Echotimeshift=0; % for multiccshot this single shot sequence does not need it
diffusion_count=size(table,1);
%% RF

tRFex=3e-3;
tRFref=5e-3;
[rf, gz, gzReph] = mr.makeSincPulse_QL_1(pi/2,'system',lims,'Duration',tRFex,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);
% Bandsep = [0 Nslices/mb]- Nslices/mb/2; % [-15, 15]
Bandsep=[0 30];
nBand   = 0;
rf_90 = mrz.makeMBPulse_less(rf.signal,'n_bands',nBand,'timeBwProduct',4,'band_sep',Bandsep);
rf.signal=rf_90;

[rf180, gz180] = mr.makeSincPulse(pi,'system',lims,'Duration',tRFref,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');
rf_180 = mrz.makeMBPulse_less(rf180.signal,'n_bands',nBand,'timeBwProduct',4,'band_sep',Bandsep);
rf180.signal=rf_180;

%% slice order

for islice_1=1:Nslices
    freqOffset_factor(islice_1)=islice_1-1-(Nslices-1)/2;
end
Nslices=Nslices/mb;
% freqOffset_factor=freqOffset_factor(Nslices/2+1:Nslices*3/2);
freqOffset_factor=freqOffset_factor(Nslices+1:Nslices*2);
%% sequence
for b_value=1:1%size(bFactor,1)
    for ds=1:5%dummy_count

        if ds ==1
            pe_enable=0;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
        else
            pe_enable=1;
        end

        for Nmulti=1:RSegment
            ro_os=2;                   % oversampling factor (in contrast to the product sequence we don't really need it)
            readoutTime=4.3e-4;        % this controls the readout bandwidth
            partFourierFactor=0.75;    % partial Fourier factor: 1: full sampling 0: start with ky=0

            % Create fat-sat pulse
            sat_ppm=-3.45;
            sat_freq=sat_ppm*1e-6*lims.B0*lims.gamma;
            rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims1,'Duration',8e-3,...
                'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
            rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
            gz_fs = mr.makeTrapezoid('z',lims1,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

            % Spoiler in front of the sequence
            spoiler_amp=3*8*42.58*10e2;
            est_rise=500e-6; % ramp time 280 us
            est_flat=2500e-6; %duration 600 us

            gp_r=mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            gp_p=mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            gp_s=mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

            gn_r=mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            gn_p=mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            gn_s=mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

            crusher_d=0.95e-3;

            gz180_crusher_1=mr.makeTrapezoid('z',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens
            gz180_crusher_2=mr.makeTrapezoid('y',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens
            gz180_crusher_3=mr.makeTrapezoid('x',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens


            gz180_c1=mr.makeTrapezoid('z',lims,'Amplitude',0*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref
            gz180_c2=mr.makeTrapezoid('y',lims,'Amplitude',-1*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref
            gz180_c3=mr.makeTrapezoid('x',lims,'Amplitude', 1*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref

            % define the output trigger to play out with every slice excitatuion
            trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

            % Define other gradients and ADC events
            deltak=1/fov;
            deltaky=RSegment*R*deltak; % Rsegement*R
            kWidth = Nx*deltak;

            % Phase blip in shortest possible time
            blip_dur = ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
            % the split code below fails if this really makes a trpezoid instead of a triangle...
            gy = mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center

            extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
            gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
            actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
            gx.amplitude=gx.amplitude/actual_area*kWidth;
            gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
            gx.flatArea = gx.amplitude*gx.flatTime;

            % calculate ADC
            adcDwellNyquist=deltak/gx.amplitude/ro_os;
            % round-down dwell time to 100 ns
            adcDwell=floor(adcDwellNyquist*1e7)*1e-7;

            adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
            % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
            adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
            % realign the ADC with respect to the gradient
            time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
            adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us

            gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims);
            [gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
            gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);

            % pe_enable support
            gy_blipup.waveform=gy_blipup.waveform*pe_enable;
            gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
            gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

            % phase encoding and partial Fourier

            Ny_pre=round((partFourierFactor-1/2)*Ny-1);  % PE steps prior to ky=0, excluding the central line
            Ny_pre=round(Ny_pre/RSegment/R);
            Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
            Ny_post=round(Ny_post/RSegment/R);
            Ny_meas=Ny_pre+Ny_post;

            % Pre-phasing gradients
            gxPre = mr.makeTrapezoid('x',lims2,'Area',-gx.area/2);
            gyPre = mr.makeTrapezoid('y',lims2,'Area',(Ny_pre*deltaky-(Nmulti-1)*R*deltak));
            [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
            % relax the PE prepahser to reduce stimulation
            gyPre = mr.makeTrapezoid('y',lims2,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
            gyPre.amplitude=gyPre.amplitude*pe_enable;

            %% gzBlip
            caipi_factor = mb;

            kz = 1/(abs(diff(Bandsep(1:2)))*thickness*caipi_factor);%
            gzBlip = mr.makeTrapezoid('z',lims,'Area',kz,'Duration',blip_dur);
            gzBlip2 = mr.makeTrapezoid('z',lims,'Area',kz*(caipi_factor-1),'Duration',blip_dur); %185 is enough

            gzBlip_parts1 = mr.splitGradientAt(gzBlip, blip_dur/2, lims);
            gzBlip_parts2 = mr.splitGradientAt(gzBlip2, blip_dur/2, lims);
            [gzBlip_l1,gzBlip_r1,~]=mr.align('right',gzBlip_parts1(1),'left',gzBlip_parts1(2),gx);
            [gzBlip_l2,gzBlip_r2,~]=mr.align('right',gzBlip_parts2(1),'left',gzBlip_parts2(2),gx);%figure,plot(gzBlip_l2.tt,gzBlip_l2.waveform)

            gzBlip_neg1=mr.addGradients({mr.scaleGrad(gzBlip_r1,-1), mr.scaleGrad(gzBlip_l1,-1)}, lims);%figure,plot(gzBlip_neg1.tt,gzBlip_neg1.waveform)
            gzBlip_pos1=mr.addGradients({mr.scaleGrad(gzBlip_r1,-1), mr.scaleGrad(gzBlip_l2,1)}, lims);%figure,plot(gzBlip_pos1.tt,gzBlip_pos1.waveform)
            gzBlip_neg2=mr.addGradients({mr.scaleGrad(gzBlip_r2,1), mr.scaleGrad(gzBlip_l1,-1)}, lims);%figure,plot(gzBlip_neg2.tt,gzBlip_neg2.waveform)


            %  gzBlipPre = mr.makeTrapezoid('z',sys,'Area',gzReph.area+kz*(caipi_factor-1)/2,'Duration',gzReph.riseTime+gzReph.flatTime+gzReph.fallTime); %mr.calcDuration(gzReph)
            gzBlipPre = mr.makeTrapezoid('z',lims,'Area', kz*(caipi_factor-1)/2,'Duration',mr.calcDuration(gyPre));
            % gzBlipReph = mr.makeTrapezoid('z',sys,'Area',gzReph.Area+kz,'Duration',gzReph.riseTime+gzReph.flatTime+gzReph.fallTime);

            gzBlipPre.amplitude=gzBlipPre.amplitude*pe_enable;
            gzBlip_neg1.waveform=gzBlip_neg1.waveform*pe_enable;gzBlip_neg1.first=gzBlip_neg1.first*pe_enable;  gzBlip_neg1.last=gzBlip_neg1.last*pe_enable;
            gzBlip_pos1.waveform=gzBlip_pos1.waveform*pe_enable;gzBlip_pos1.first=gzBlip_pos1.first*pe_enable;  gzBlip_pos1.last=gzBlip_pos1.last*pe_enable;
            gzBlip_neg2.waveform=gzBlip_neg2.waveform*pe_enable;gzBlip_neg2.first=gzBlip_neg2.first*pe_enable;  gzBlip_neg2.last=gzBlip_neg2.last*pe_enable;
            gzBlip_l1.waveform=gzBlip_l1.waveform*pe_enable; gzBlip_l1.first=gzBlip_l1.first*pe_enable;  gzBlip_l1.last=gzBlip_l1.last*pe_enable;
            gzBlip_l2.waveform=gzBlip_l2.waveform*pe_enable;gzBlip_l2.first=gzBlip_l2.first*pe_enable;  gzBlip_l2.last=gzBlip_l2.last*pe_enable;
            gzBlip_r1.waveform=gzBlip_r1.waveform*pe_enable;gzBlip_r1.first=gzBlip_r1.first*pe_enable;  gzBlip_r1.last=gzBlip_r1.last*pe_enable;
            gzBlip_r2.waveform=gzBlip_r2.waveform*pe_enable;gzBlip_r2.first=gzBlip_r2.first*pe_enable;  gzBlip_r2.last=gzBlip_r2.last*pe_enable;

            %% for after the first echo QL
            gxPre_post = mr.makeTrapezoid('x',lims2,'Area',-gx.area/2);
            gyPre_post = mr.makeTrapezoid('y',lims2,'Area',(Ny_post-1)*deltaky); % QL

            [gxPre_post,gyPre_post]=mr.align('right',gxPre_post,'left',gyPre_post);
            % relax the PE prepahser to reduce stimulation
            gyPre_post = mr.makeTrapezoid('y',lims2,'Area',gyPre_post.area,'Duration',mr.calcDuration(gxPre_post,gyPre_post));

            gyPre_post.amplitude=gyPre_post.amplitude*pe_enable;

            if (mod(Ny_meas,2)==0)
                gxPre_post.amplitude=-gxPre_post.amplitude;
            end

            gzBlip_post=mr.makeTrapezoid('z',lims2,'Area', -kz*(caipi_factor-1)/2,'Duration',mr.calcDuration(gxPre_post));
            gzBlip_post.amplitude=gzBlip_post.amplitude*pe_enable;

            if (ds==3)
                gyPre=mr.scaleGrad(gyPre,-1);
                gy_blipdownup=mr.scaleGrad(gy_blipdownup,-1);
                gy_blipdown=mr.scaleGrad(gy_blipdown,-1);
                gy_blipup=mr.scaleGrad(gy_blipup,-1);
                gyPre_post=mr.scaleGrad(gyPre_post,-1);
            end

            if (ds < 4) % need to adjust delta TE1 and delta TE2 for b0, bcuz of the crusher
                crusher_switch=1;
            else
                crusher_switch=0;
            end

            % Calculate delay times
            durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
            rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
            rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);

            delayTE1=ceil((TE1/2 - mr.calcDuration(rf,gz) - mr.calcDuration(gzReph) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
            delayTE1=delayTE1-crusher_switch*mr.calcDuration(gz180_crusher_1); %QL
            delayTE2tmp=ceil((TE1/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
            assert(delayTE1>=0);

            gxPre.delay=0;
            gyPre.delay=0;
            delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre);
            delayTE2=delayTE2-crusher_switch*mr.calcDuration(gz180_crusher_1); %QL
            [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
            assert(delayTE2>=0);

            % match tdm
            dur_match=0.01438;
            delayTE2=delayTE2-dur_match;

            dmatch=mr.makeDelay(dur_match);
            dur_match_TE_1=0.5*(TE2-TE1);
            dmatch_TE_1=mr.makeDelay(dur_match_TE_1);
            dur_match_TE_2=0.5*(TE3-TE1);
            dmatch_TE_2=mr.makeDelay(dur_match_TE_2);

            dur2end=(Ny_post-0.5)*mr.calcDuration(gx)+mr.calcDuration(gxPre_post,gyPre_post);

            small_delta=delayTE2-ceil(lims2.maxGrad/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            big_delta=delayTE1+mr.calcDuration(rf180,gz180)+0.0027+2*mr.calcDuration(rf180,gz180);% to match TDM

            g=sqrt(bFactor(b_value,ds)*1e6/bFactCalc(1,small_delta,big_delta)); % QL
            gr=ceil(g/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims2);
            %% Qiang: split the diffusion gradient to 3 axis, new version based on Jon and Maxim's suggestions

            g_x=g.*table(ds,1,b_value);
            g_y=g.*table(ds,2,b_value);
            g_z=g.*table(ds,3,b_value);

            if (table(ds,1,b_value )==0&&table(ds,2, b_value)==0&&table(ds,3, b_value)==0) % when the b-value=0, the section using Sepherical coorodinate system is wrong.
                gDiff_x=gDiff; gDiff_x.channel='x';
                gDiff_y=gDiff; gDiff_y.channel='y';
                gDiff_z=gDiff; gDiff_z.channel='z';

            else
                [azimuth,elevation,r] = cart2sph(g_x,g_y,g_z);
                polar= -(pi/2-elevation);

                Gr=mr.rotate('z',azimuth,mr.rotate('y',polar,gDiff));
                if size(Gr,2)==3
                    gDiff_x=Gr{1,2};
                    gDiff_y=Gr{1,3};
                    gDiff_z=Gr{1,1};
                else
                    if size(Gr,2)==2
                        diffusion_blank=find( table(ds,:, b_value)==0);
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
            end


            duration_diff=mr.calcDuration(gDiff);

            assert(duration_diff<=delayTE1);
            assert(duration_diff<=delayTE2);

            %% Calculate the echo time shift for multishot EPI Qiang Liu
            actual_esp=gx.riseTime+gx.flatTime+gx.fallTime;

            %% TR calculation
            TR=2500e-3;
            delayTR=TR/Nslices*3-mr.calcDuration(gp_r) - mr.calcDuration(gn_r) - TE - rfCenterInclDelay-dur2end; % *2
            delayTR=round(delayTR,3);
            dTR=mr.makeDelay(delayTR);

            %% Define sequence blocks
            for s=1:3:Nslices
                seq.addBlock(gp_r,gp_p,gp_s);
                seq.addBlock(rf_fs, gn_r,gn_p,gn_s);

                rf.freqOffset=gz.amplitude*thickness*freqOffset_factor(s);
                %     rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
                rf180.freqOffset=gz180.amplitude*thickness*freqOffset_factor(s);
                %     rf180.phaseOffset=pi/2-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
                rf180.phaseOffset=pi/2;
                seq.addBlock(rf,gz,trig);
                seq.addBlock(gzReph);

                switch TE_index
                    case 1
                        %
                    case 2
                        seq.addBlock(dmatch_TE_1);
                    case 3
                        seq.addBlock(dmatch_TE_2);
                end
                seq.addBlock(mr.makeDelay(delayTE1),gDiff_x, gDiff_y, gDiff_z);

                if (ds < 4)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end

                seq.addBlock(rf180,gz180); % QL: it used to be gz180n

                if (ds < 4)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
                seq.addBlock(dmatch);% to match tdm
                seq.addBlock(mr.makeDelay(delayTE2),gDiff_x, gDiff_y, gDiff_z);

                switch TE_index
                    case 1
                        %
                    case 2
                        seq.addBlock(dmatch_TE_1);
                    case 3
                        seq.addBlock(dmatch_TE_2);
                end

                if ds==1
                    gy_blipup.last=0;
                    gy_blipdown.first=0;
                    gy_blipdownup.last=0;
                    gy_blipdownup.first=0;
                end

                if cai_flag
                    gzPre = gzBlipPre;%gzReph
                end
                seq.addBlock(mr.align('left',gyPre,gzPre,'right',gxPre));
                %% Blip-CAIPI
                for i=1:Ny_meas % loop over phase encodes
                    gx_temp = mr.scaleGrad(gx,(-1)^(i-1));
                    if i==1
                        if cai_flag
                            step = 0;
                            seq.addBlock(gx_temp,gy_blipup,mr.scaleGrad(gzBlip_l1,-1),adc); % Read the first line of k-space with a single half-blip at the end
                        else
                            seq.addBlock(gx_temp,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                        end
                    elseif i==Ny_meas
                        if cai_flag
                            step = step + 1;
                            %check the gz moment
                            switch mod(step,caipi_factor)
                                case num2cell(1:caipi_factor-2)
                                    seq.addBlock(gx_temp,gy_blipdown,mr.scaleGrad(gzBlip_r1,-1),adc);
                                case caipi_factor-1
                                    seq.addBlock(gx_temp,gy_blipdown,mr.scaleGrad(gzBlip_r1,-1),adc);
                                case 0
                                    seq.addBlock(gx_temp,gy_blipdown,mr.scaleGrad(gzBlip_r2,1),adc);
                            end
                        else
                            seq.addBlock(gx_temp,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                        end
                    else
                        if cai_flag
                            step = step+1;
                            switch mod(step,caipi_factor)
                                case num2cell(1:caipi_factor-2)
                                    seq.addBlock(gx_temp,gy_blipdownup,gzBlip_neg1,adc);
                                case caipi_factor-1
                                    seq.addBlock(gx_temp,gy_blipdownup,gzBlip_pos1,adc);
                                case 0
                                    seq.addBlock(gx_temp,gy_blipdownup,gzBlip_neg2,adc);
                            end
                        else
                            seq.addBlock(gx_temp,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    end
                end
                seq.addBlock(gxPre_post,gyPre_post,gzBlip_post);
                seq.addBlock(dTR); % seperate
            end %slice loop
        end % for multishot
        disp(['diffusion dir: ', num2str(ds)])
    end % EPI REF diffusion LOOP
end% b-value
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

% seq.plot();             % Plot sequence waveforms

% trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
%
% % plot k-spaces
% figure; plot(t_ktraj, ktraj'); % plot cthe entire k-space trajectory
% hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi-diff');
%  [pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/MP_GPA_K2309_2250V_951A_AS82.asc'); % TERRA-XR

%  if pns_ok

seq.write(seq_file);
%  end

% seq.install('siemens');
% seq.sound(); % simulate the seq's tone

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits
%
% rep = seq.testReport;
% fprintf([rep{:}]);
%
%%
function b=bFactCalc(g, delta, DELTA)
sigma=1;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end
