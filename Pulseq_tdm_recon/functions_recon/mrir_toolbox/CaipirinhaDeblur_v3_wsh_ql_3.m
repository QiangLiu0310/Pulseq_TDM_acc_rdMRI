function [Kcorrected] = CaipirinhaDeblur_v3_wsh_ql_3(K)
% 
%  load('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/data_processing/data_processing_01_07/slicepos.mat')
%  PhaseShiftPerMM = (pi/65); 
% load('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/data_processing/data_processing_01_14/slicepos_offcenter.mat')
% PhaseShiftPerMM = (pi/13); % pulseq
% PhaseShiftPerMM = (2*pi/3/20); % pulseq
PhaseShiftPerMM = (2*pi/2/6); % pulseq

% SlicePos=[-12.5:1:-0.5]; % pulseq


 SlicePos=[5.5:-1:0];

% SlicePos  = SlicePos(end:-1:1); 

%  SlicePos=[-29.5:1:29.5];
PhaseShiftBtwSimulSlices= 2*pi/2 ;%2
Kcorrected = (zeros(size(K)));

% if PhaseShiftBtwSimulSlices ~= 0
    for SlcCount = 1: size(K,3), % ( length(prot.sSliceArray.asSlice) > 1 ),
      PhaseShift = PhaseShiftPerMM*SlicePos(SlcCount);
      % fprintf(' %f',PhaseShift);
        if abs(PhaseShiftBtwSimulSlices) == pi % FOV/2 shift
            Kcorrected(:,2:2:end,SlcCount,:) =  K(:,2:2:end,SlcCount,:);
            Kcorrected(:,1:2:end,SlcCount,:) =  K(:,1:2:end,SlcCount,:)*exp(-i*PhaseShift);
         elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/3 % FOV/3 shift
            Kcorrected(:,1:3:end,SlcCount,:) =  K(:,1:3:end,SlcCount,:);
            Kcorrected(:,2:3:end,SlcCount,:) =  K(:,2:3:end,SlcCount,:)*exp(-i*PhaseShift);
            Kcorrected(:,3:3:end,SlcCount,:) =  K(:,3:3:end,SlcCount,:)*exp(-i*2*PhaseShift);

        elseif abs(PhaseShiftBtwSimulSlices) == pi/2 % FOV/4 shift
            Kcorrected(:,1:4:end,SlcCount,:) =  K(:,1:4:end,SlcCount,:);
            Kcorrected(:,2:4:end,SlcCount,:) =  K(:,2:4:end,SlcCount,:)*exp(-i*PhaseShift);
            Kcorrected(:,3:4:end,SlcCount,:) =  K(:,3:4:end,SlcCount,:)*exp(-i*2*PhaseShift);
            Kcorrected(:,4:4:end,SlcCount,:) =  K(:,4:4:end,SlcCount,:)*exp(-i*3*PhaseShift);
        end
%     end
    % fprintf('\n');
% else
%     Kcorrected = K;
    end
% end



