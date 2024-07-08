%% based on epi diffusion rs @Pulseq | @Siemens code by Yang Ji Phd
%% TDM-EPI 3 echo different echo times
%% GRAPPA=3
%% TE= 72/102/132 ms
%%  introduced components
%    echo time shift gradients
%    the ss gradient for 2nd 180 pulse was reversed
%% for magnus
%    Diffusion slew rate 176.6 T/m/s grad mag: 284.776 mT/m
%    add labels (April 17 2024)
% reduce the slew rate for the gy post
%% Qiang Liu Jan/17/2023 at PNL

clear all;close all;clc
field=3.00; % prisma 2.89; magnus 3.00
seq_name='tdm_3e_mb1_21sli_2mm_R3_test.seq';

lims = mr.opts('MaxGrad',50,'GradUnit','mT/m',...
    'MaxSlew',740,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', field);
lims1 = mr.opts('MaxGrad',100,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', field);
lims2 = mr.opts('MaxGrad',284,'GradUnit','mT/m',...
    'MaxSlew',176,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', field);

%%  GE  system QL
lims.rfRingdownTime=6*1e-05; % 60 us
lims1.rfRingdownTime=6*1e-05;
lims1.adcDeadTime=20e-6; % EPI
lims.adcDeadTime=20e-6;
commonRasterTime = 20e-6;
%% Pulseq sequence

tmp=load('/Users/ln915/Documents/code/Pulseq_Qiang_PNL/matlab/demoSeq/diffusion_table/TMD_part1_bvec.txt');
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
fov=220e-3; Nx=110; Ny=Nx; % Define FOV and resolution
thickness=2.0e-3;            % slice thinckness
Nslices=21;
TE1=55e-3;
TE2=80e-3;
TE3=105e-3;
TE_incre=6.0e-3; % the second set of TE: [82 107 132]
tdm=3;
R=3;
est_rise=280e-6; % ramp time 280 us
est_flat=800e-6; %duration 600 us
amp_factor=20*42.58*10e2; % kHz/m for crusher
%% interleaved slices
for islice_1=1:Nslices
    freqOffset_factor(islice_1)=islice_1-1-(Nslices-1)/2;
end
% Nslices=Nslices/mb;
% freqOffset_factor=freqOffset_factor(Nslices/2+1:Nslices*3/2);
for b_value=1:size(bFactor,1)
    for tdmloop=1: tdm

        slice_index=[1:Nslices]; slice_index=reshape(slice_index,[Nslices/tdm, tdm])';
        slice_index1=zeros(tdm,Nslices/tdm);

        switch tdmloop
            case 1
                slice_index=reshape(slice_index,[Nslices,1]);
                interleaved_freqOffset_factor=freqOffset_factor(slice_index);
            case 2
                slice_index1(1,:)=slice_index(2,:); slice_index1(2,:)=slice_index(3,:); slice_index1(3,:)=slice_index(1,:);
                slice_index=reshape(slice_index1,[Nslices,1]);
                interleaved_freqOffset_factor=freqOffset_factor(slice_index);
            case 3
                slice_index1(1,:)=slice_index(3,:); slice_index1(2,:)=slice_index(1,:); slice_index1(3,:)=slice_index(2,:);
                slice_index=reshape(slice_index1,[Nslices,1]);
                interleaved_freqOffset_factor=freqOffset_factor(slice_index);
        end
        slice_index=reshape(slice_index,[Nslices,1]);
        interleaved_freqOffset_factor=freqOffset_factor(slice_index);

        for ds=1:3%dummy_count % actually it is diffusion count

            if ds == 1
                pe_enable=0;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
            else
                pe_enable=1;
            end              % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
            ro_os=1;                   % oversampling factor (in contrast to the product sequence we don't really need it)
            readoutTime=2.8e-4;        % this controls the readout bandwidth  2000/0.58 ms % 3.9
            partFourierFactor=6/8;    % partial Fourier factor: 1: full sampling 0: start with ky=0

            tRFex=3e-3;
            tRFref=5e-3;

            % Create fat-sat pulse
            sat_ppm=-3.45;
            sat_freq=sat_ppm*1e-6*lims.B0*lims.gamma;
            rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims1,'Duration',8e-3,...
                'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
            rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
            gz_fs = mr.makeTrapezoid('z',lims1,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

            % Spoiler in front of the sequence
            spoiler_amp=3*8*42.58*10e2;
            est_rise_p=500e-6; % ramp time 280 us
            est_flat_p=2500e-6; %duration 600 us

            gp_r=mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise_p,'flatTime',est_flat_p,'system',lims1);
            gp_p=mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise_p,'flatTime',est_flat_p,'system',lims1);
            gp_s=mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise_p,'flatTime',est_flat_p,'system',lims1);

            gn_r=mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise_p,'flatTime',est_flat_p,'system',lims1);
            gn_p=mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise_p,'flatTime',est_flat_p,'system',lims1);
            gn_s=mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise_p,'flatTime',est_flat_p,'system',lims1);



            % Create 90 degree slice selection pulse and gradient
            [rf1, gz1, gzReph,gzReph_double] = mr.makeSincPulse_QL_1_tdm(pi/2,'system',lims,'Duration',tRFex,...
                'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);
            rf2=rf1;gz2=gz1;
            rf3=rf1;gz3=gz1;

            % Create 180 degree slice refocusing pulse and gradients
            [rf180, gz180] = mr.makeSincPulse(pi,'system',lims,'Duration',tRFref,...
                'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');
            % gz180_r=mr.makeTrapezoid('z',lims,'Area',-gz180.area,'Duration',est_rise*2+est_flat);
            gz180_r=mr.scaleGrad(gz180,-1);
            rf2180=rf180;
            rf3180=rf180;


            crusher_d=0.95e-3;
            gz180_crusher_1=mr.makeTrapezoid('z',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens
            gz180_crusher_2=mr.makeTrapezoid('y',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens
            gz180_crusher_3=mr.makeTrapezoid('x',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens


            % define the output trigger to play out with every slice excitatuion
            trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

            % Define other gradients and ADC events
            deltak=1/fov;
            deltaky=R*deltak;
            kWidth = Nx*deltak;

            % Phase blip in shortest possible time
            blip_dur = ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
            gy = mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center

            % add for GE QL
            gy=trap4ge(gy,commonRasterTime,lims);

            extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
            gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
            actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
            gx.amplitude=gx.amplitude/actual_area*kWidth;
            gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
            gx.flatArea = gx.amplitude*gx.flatTime;

            % calculate ADC
            adcDwellNyquist=deltak/gx.amplitude/ro_os;
            %             adcDwellNyquist=2e-06;
            % round-down dwell time to 100 ns
            adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
            adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
            % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
            adc = mr.makeAdc(adcSamples,'system', lims, 'Dwell',adcDwell,'Delay',blip_dur/2);
            % realign the ADC with respect to the gradient
            time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
            adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us

            % split the blip into two halves and produce a combined synthetic gradient
            gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims);
            [gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
            gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);

            % pe_enable support
            gy_blipup.waveform=gy_blipup.waveform*pe_enable;
            gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
            gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

            % phase encoding and partial Fourier
            Ny_pre=round((partFourierFactor-1/2)*Ny-1);  % PE steps prior to ky=0, excluding the central line
            Ny_pre=round(Ny_pre/R);
            Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
            Ny_post=round(Ny_post/R);
            Ny_meas=Ny_pre+Ny_post;

            % Pre-phasing gradients
            gxPre = mr.makeTrapezoid('x',lims2,'Area',-gx.area/2);
            gyPre = mr.makeTrapezoid('y',lims2,'Area',Ny_pre*deltaky);
            [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
            % relax the PE prepahser to reduce stimulation
            gyPre = mr.makeTrapezoid('y',lims2,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
            gyPre.amplitude=gyPre.amplitude*pe_enable;

            gxPre_post = mr.makeTrapezoid('x',lims2,'Area',-gx.area/2);
            gyPre_post = mr.makeTrapezoid('y',lims2,'Area',(Ny_post-1)*deltaky);

            [gxPre_post,gyPre_post]=mr.align('right',gxPre_post,'left',gyPre_post);
            % relax the PE prepahser to reduce stimulation
            gyPre_post = mr.makeTrapezoid('y',lims2,'Area',gyPre_post.area,'Duration',mr.calcDuration(gxPre_post,gyPre_post));
            gyPre_post.amplitude=gyPre_post.amplitude*pe_enable;

            if (mod(Ny_meas,2)==0)
                gxPre_post.amplitude=-gxPre_post.amplitude;
            end
            % % add for GE QL
            % gyPre_post=trap4ge(gyPre_post,commonRasterTime,lims2);
            % gxPre_post=trap4ge(gxPre_post,commonRasterTime,lims2);
            % duration of EPI mod
            epi_dur=mr.calcDuration(gx)*Ny_meas+mr.calcDuration(gxPre_post,gyPre_post)+mr.calcDuration(gyPre);

            if (ds==3)
                gyPre=mr.scaleGrad(gyPre,-1);
                gy_blipdownup=mr.scaleGrad(gy_blipdownup,-1);
                gy_blipdown=mr.scaleGrad(gy_blipdown,-1);
                gy_blipup=mr.scaleGrad(gy_blipup,-1);
                gyPre_post=mr.scaleGrad(gyPre_post,-1);
            end
            % add for ge
            gyPre=trap4ge(gyPre,commonRasterTime,lims);
            gxPre=trap4ge(gxPre,commonRasterTime,lims);

            if (ds < 4) % need to adjust delta TE1 and delta TE2 for b0, bcuz of the crusher
                crusher_switch=1;
            else
                crusher_switch=0;
            end

            % Calculate delay times
            d2_s_dur=est_rise*2+est_flat;
            dm1_s_dur=d2_s_dur;
            delay2=0.5*(TE2-TE1)-epi_dur+mr.calcDuration(rf180,gz180);
            d_delay2=mr.makeDelay(delay2);
            delay3=0.5*(TE3-TE2)-epi_dur+mr.calcDuration(rf180,gz180);
            d_delay3=mr.makeDelay(delay3);

            durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
            rfCenterInclDelay=rf1.delay + mr.calcRfCenter(rf1);
            rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
            delayTE1=ceil((TE1/2 - mr.calcDuration(rf1,gz1) + rfCenterInclDelay- mr.calcDuration(gzReph) - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
            delayTE1=delayTE1-crusher_switch*mr.calcDuration(gz180_crusher_1);
            %         delayTE1=delayTE1-dur_rephase; % not useful

            delayTE2tmp=ceil((TE1/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
            assert(delayTE1>=0);
            gxPre.delay=0;gyPre.delay=0;
            delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre);
            delayTE2=delayTE2-crusher_switch*mr.calcDuration(gz180_crusher_1);
            [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);

            % delayTE2=delayTE2-d2_s_dur-dm1_s_dur-epi_dur-rf_distance; % tdm
            % delayTE2=delayTE2-d2_s_dur-dm1_s_dur-epi_dur; % tdm
            delayTE2=delayTE2-d2_s_dur*3-mr.calcDuration(rf2180,gz180)-mr.calcDuration(rf2180,gz180);% tdm

            assert(delayTE2>=0);

            small_delta=delayTE2-ceil(lims2.maxGrad/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            % big_delta=delayTE1+mr.calcDuration(rf180,gz180);
            % big_delta=delayTE1+mr.calcDuration(rf180,gz180)+mr.calcDuration(rf180,gz180)+d2_s_dur+rf_distance_180_2;
            big_delta=delayTE1+mr.calcDuration(rf180,gz180)+2*d2_s_dur+2*mr.calcDuration(rf180,gz180);
            g=sqrt(bFactor(b_value,ds)*1e6/bFactCalc(1,small_delta,big_delta));
            gr=ceil(g/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;

            gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims2); 
            % outside the sequence loop 
            % inside: only scaling (b-factor)
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

            %% crushers for echo time shift
            crudi_r=0;crudi_p=0;crudi_s=0;

            if (g_y==0)
                crudi_s=0;
                crudi_p=1;
                crudi_r=0;
            else
                crudi_s=0;
                crudi_p=1/sqrt(1+(g_y^2)./(g_x)^2);
                crudi_r=-g_y/g_x*crudi_p;
            end

            g1_r=mr.makeTrapezoid('x','amplitude',crudi_r*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            g1_p=mr.makeTrapezoid('y','amplitude',crudi_p*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            % gm1_r=mr.makeTrapezoid('x','amplitude',-crudi_r*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims);
            % gm1_p=mr.makeTrapezoid('y','amplitude',-crudi_p*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims);
            g2_r=mr.makeTrapezoid('x','amplitude',2*crudi_r*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            g2_p=mr.makeTrapezoid('y','amplitude',2*crudi_p*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            gm2_r=mr.makeTrapezoid('x','amplitude',-2*crudi_r*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
            gm2_p=mr.makeTrapezoid('y','amplitude',-2*crudi_p*amp_factor,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

            %% calculate the delay for 90 ex
            delay1=0.5*(TE2-TE1)-d2_s_dur-mr.calcDuration(rf180,gz180)-mr.calcDuration(rf1,gz1)-mr.calcDuration(gzReph_double);
            d_delay1=mr.makeDelay(delay1);

            %% TR calculation Qiang Liu
            TR=2500e-3;
            Nsli=Nslices/tdm;
            delayTR=ceil((TR-Nsli*mr.calcDuration(gp_r) - Nsli*mr.calcDuration(gn_r)...
                -TE3*Nsli-Nsli*Ny_post*mr.calcDuration(gx)-Nsli*mr.calcDuration(gxPre_post)-Nsli*mr.calcDuration(rf1,gz1)*0.5...
                )/seq.gradRasterTime)*seq.gradRasterTime;
            delayTR=delayTR/Nsli;
            delayTR=round(delayTR,3); % pulseq will allow delay time at power -3
            % delayTR=0.01;
            dTR=mr.makeDelay(delayTR);
            % segmentID=1;% only add 1 segmentID to reference scan
            %% Define sequence blocks
            for s=1:tdm:3%Nslices

                % if (s==1)&&(ds==1)&&(tdmloop==1)&&(b_value==1)
                seq.addBlock(gp_r,gp_p,gp_s,mr.makeLabel('SET', 'TRID', ds)); % TRID=ds
                % else
                %     seq.addBlock(gp_r,gp_p,gp_s);
                % end
                seq.addBlock(rf_fs, gn_r,gn_p,gn_s);

                rf1.freqOffset=gz1.amplitude*thickness*interleaved_freqOffset_factor(s);
                rf180.freqOffset=gz180.amplitude*thickness*interleaved_freqOffset_factor(s);
                rf180.phaseOffset=pi/2;

                rf2.freqOffset=gz2.amplitude*thickness*interleaved_freqOffset_factor(s+1);
                rf2180.freqOffset=-gz180.amplitude*thickness*interleaved_freqOffset_factor(s+1); % -
                rf2180.phaseOffset=pi/2;

                rf3.freqOffset=gz3.amplitude*thickness*interleaved_freqOffset_factor(s+2);
                %                 rf3180.freqOffset=gz180.amplitude*thickness*interleaved_freqOffset_factor(s+2);
                rf3180.phaseOffset=pi/2;

                seq.addBlock(rf3,gz3);
                seq.addBlock(gzReph_double);
                seq.addBlock(d_delay1);
                seq.addBlock(rf2,gz2);
                seq.addBlock(gzReph_double); %rewinder QL
                seq.addBlock(d_delay1);
                seq.addBlock(rf1,gz1);
                seq.addBlock(gzReph);
                %             seq.addBlock(dphase); % SPlit the rewinder and the diffusion gradient
                seq.addBlock(mr.makeDelay(delayTE1),gDiff_x, gDiff_y, gDiff_z);

                if (ds < 4)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
                seq.addBlock(rf180,gz180);
                seq.addBlock(g1_r,g1_p);
                seq.addBlock(rf2180,gz180_r);
                seq.addBlock(g1_r,g1_p);
                seq.addBlock(rf3180,gz180);

                if (ds < 4)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
                seq.addBlock(mr.makeDelay(delayTE2),gDiff_x, gDiff_y, gDiff_z);

                seq.addBlock(gm2_r,gm2_p);

                if ds ==1
                    gy_blipup.last=0;
                    gy_blipdown.first=0;
                    gy_blipdownup.last=0;
                    gy_blipdownup.first=0;
                end

                %EPI For 3rd slice
                seq.addBlock(gxPre,gyPre);
                for i=1:Ny_meas
                    if i==1
                        seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                    elseif i==Ny_meas
                        seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                    else
                        seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    end
                    %         gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    if (mod(Ny_meas,2)==1)
                        if (i<Ny_meas)
                            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                        end
                    else
                        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    end
                end
                seq.addBlock(gxPre_post,gyPre_post);

                seq.addBlock(g2_r,g2_p); % ETS gradient
                seq.addBlock(d_delay2);
                %EPI For 2nd slice
                seq.addBlock(gxPre,gyPre);
                for i=1:Ny_meas
                    if i==1
                        seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                    elseif i==Ny_meas
                        seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                    else
                        seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    end
                    %         gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    if (mod(Ny_meas,2)==1)
                        if (i<Ny_meas)
                            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                        end
                    else
                        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    end
                end

                seq.addBlock(gxPre_post,gyPre_post);

                seq.addBlock(g2_r,g2_p); % ETS gradient
                seq.addBlock(d_delay3);
                % EPI for 1st slice
                seq.addBlock(gxPre,gyPre);
                for i=1:Ny_meas
                    if i==1
                        seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                    elseif i==Ny_meas
                        seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                    else
                        seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    end
                    %         gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    if (mod(Ny_meas,2)==1)
                        if (i<Ny_meas)
                            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                        end
                    else
                        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    end
                end
                seq.addBlock(gxPre_post,gyPre_post);
                seq.addBlock(dTR);
            end%slice
        end %REF
    end % TDM
end
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
% plot k-spaces
% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi-diff');
% [pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/MP_GPA_K2309_2250V_951A_AS82.asc'); % TERRA-XR

seq.write(seq_name); % v2 reversed the freqoffset of the second 180 pulse, just for test
%%
function b=bFactCalc(g, delta, DELTA)
sigma=1;

kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end
