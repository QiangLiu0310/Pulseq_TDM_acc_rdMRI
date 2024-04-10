function [Kcorrected] = CaipirinhaDeblur_v3_wsh(K, prot, evp , PhaseShiftBtwSimulSlices, SliceSep)
%
% Kcorrected = CaipirinhaDeblur_v3_wsh(K, prot, evp, [PhaseShift, SliceSep]);
%
% K    : sampled SMS kspace
% prot : protocol structure, populated by eval_ascconv
% evp  : can set evp.NSlcMeas to number of uncollapsed slices

% optional (if unspecified, will be extracted from prot)
%
% PhaseShift : the FOV shift applied to the collapsed slices
% SliceSep   : the physical distance (in mm) between the collapsed slices

% Kawin Setsompop 
% Oct 19 2009

% now work for all acquisition orientations
% assume that the slices are already sorted (so can use para in original prot to calc stuff) 

% SliceSep in mm

if (nargin < 4),
  [~,~,~,~, PhaseShiftBtwSimulSlices, SliceSep] = mrir_array_SMS_recon_params(prot,evp);
end;

% correction of SlicePos value due to rotation of principle gradient axis
if size(K,3) > 1;

    % When the value for these coordinate is very small, e.g.
    % sPosition_dSag = 1.343452e-9 then the read_meas_dat function that readin the data would not
    % recognize it and will leave the array empty so fix it here
    if ~isfield(prot.sSliceArray.asSlice(1).sPosition,'dSag') && ...
          ~isfield(prot.sSliceArray.asSlice(2).sPosition,'dSag')
        for count = 1:length(prot.sSliceArray.asSlice)
            prot.sSliceArray.asSlice(count).sPosition.dSag = 0;
        end
    end
    if ~isfield(prot.sSliceArray.asSlice(1).sNormal,'dSag') && ...
          ~isfield(prot.sSliceArray.asSlice(2).sNormal,'dSag')
        for count = 1:length(prot.sSliceArray.asSlice)
            prot.sSliceArray.asSlice(count).sNormal.dSag = 0;
        end
    end
    if ~isfield(prot.sSliceArray.asSlice(1).sPosition,'dCor') && ...
          ~isfield(prot.sSliceArray.asSlice(2).sPosition,'dCor')
      %  isempty(prot.sSliceArray(1).sPosition_dCor) && isempty(prot.sSliceArray(2).sPosition_dCor)
        for count = 1:length(prot.sSliceArray.asSlice)
            prot.sSliceArray.asSlice(count).sPosition.dCor = 0;
        end
    end
    if ~isfield(prot.sSliceArray.asSlice(1).sNormal,'dCor') && ...
          ~isfield(prot.sSliceArray.asSlice(2).sNormal,'dCor')
      %  isempty(prot.sSliceArray(1).sPosition_dCor) && isempty(prot.sSliceArray(2).sPosition_dCor)
        for count = 1:length(prot.sSliceArray.asSlice)
            prot.sSliceArray.asSlice(count).sNormal.dCor = 0;
        end
    end
    if ~isfield(prot.sSliceArray.asSlice(1).sPosition,'dTra') && ...
          ~isfield(prot.sSliceArray.asSlice(2).sPosition,'dTra')
      % isempty(prot.sSliceArray(1).sPosition_dTra) && isempty(prot.sSliceArray(2).sPosition_dTra)
        for count = 1:length(prot.sSliceArray.asSlice)
            prot.sSliceArray.asSlice(count).sPosition.dTra = 0;
        end
    end
    if ~isfield(prot.sSliceArray.asSlice(1).sNormal,'dTra') && ...
          ~isfield(prot.sSliceArray.asSlice(2).sNormal,'dTra')
        for count = 1:length(prot.sSliceArray.asSlice)
            prot.sSliceArray.asSlice(count).sNormal.dTra = 0;
        end
    end
    
    NormalVec = [ prot.sSliceArray.asSlice(1).sNormal.dSag, ...
                  prot.sSliceArray.asSlice(1).sNormal.dCor, ...
                  prot.sSliceArray.asSlice(1).sNormal.dTra].';   
    for cnt = 1 : length(prot.sSliceArray.asSlice),
      Pos(cnt,1) = [prot.sSliceArray.asSlice(cnt).sPosition.dSag].';
      Pos(cnt,2) = [prot.sSliceArray.asSlice(cnt).sPosition.dCor].';
      Pos(cnt,3) = [prot.sSliceArray.asSlice(cnt).sPosition.dTra].';
    end;
    SlicePos = Pos*NormalVec;
else
    keyboard('only single slice data so cant determine correctionFac if Gz is rotated')
end
      
%SlicePos  = SlicePos(end:-1:1); 
PhaseShiftPerMM = (PhaseShiftBtwSimulSlices/SliceSep);

%Kcorrected = single(zeros(size(K)));
Kcorrected = (zeros(size(K)));

if PhaseShiftBtwSimulSlices ~= 0
    for SlcCount = 1: size(K,3), % ( length(prot.sSliceArray.asSlice) > 1 ),
      PhaseShift = PhaseShiftPerMM*SlicePos(SlcCount);
      % fprintf(' %f',PhaseShift);
        if abs(PhaseShiftBtwSimulSlices) == pi % FOV/2 shift
            Kcorrected(:,1:2:end,SlcCount,:) =  K(:,1:2:end,SlcCount,:);
            Kcorrected(:,2:2:end,SlcCount,:) =  K(:,2:2:end,SlcCount,:)*exp(-i*PhaseShift);
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
    end
    % fprintf('\n');
else
    Kcorrected = K;
end



