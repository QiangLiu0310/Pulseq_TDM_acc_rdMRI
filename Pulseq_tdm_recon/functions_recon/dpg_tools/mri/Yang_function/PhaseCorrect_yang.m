function [Raw_image_pc, phaseCoeff] = PhaseCorrect_yang(Raw_PC,Raw_image)

% % Raw_PC dimension: Nfreq  Ncoil  ROp/ROn  Ntotal=Nslice*Nbvalue*Ndirection
% [Nfreq,  Ncoil, Necho, Ntotal]=size(Raw_PC);
% 
% indx = 1:(Nfreq-1);
% 
% ROp_data=squeeze(Raw_PC(:,:,1,:));
% ROn_data=squeeze(Raw_PC(:,:,2,:));
% 
% % the simplest model of NGC employs a linear + scalar phase offset, which
% % manifests as shift in the image domain.  Convert to (hybrid) image domain here:
% ROp=ffts(ifft(ffts(ROp_data,1),[],1),1);
% ROn=ffts(ifft(ffts(ROn_data,1),[],1),1);
% 
% % compute the NGC correction parameters (for all coils and Ntotal(Nslice*Nbvalue*Ndirect)
% P0 = ROn(indx+0,:).*conj(ROp(indx+0,:));
% P1 = ROn(indx+1,:).*conj(ROp(indx+1,:));
% figure;plot(angle(P0))
% 
% epi_pc_lin_tem = angle( sum( P1.* conj(P0), 1 ));
% epi_pc_lin = reshape(epi_pc_lin_tem,[Ncoil Ntotal]);
% phaseCoeff(2,:,:) = epi_pc_lin;
% 
% % apply the linear phase term to the phase difference signals,
% x_index = ([1:Nfreq]-Nfreq/2)';
% C0t = x_index*epi_pc_lin_tem;
% C0tt=reshape(C0t,[Nfreq Ncoil Ntotal]);
% C0ttt = ROp.* conj(ROn) .* exp(j*C0tt);
% temp=C0ttt(:,:);
% figure;plot(angle(temp));
% 
% epi_pc_con = angle(sum(C0ttt,1));
% phaseCoeff(1,:,:) = epi_pc_con;
%% scott
% Raw_PC dimension: Nfreq  Ncoil  ROp/ROn  Ntotal=Nslice*Nbvalue*Ndirection
[Nfreq,  Ncoil, NPhaCorrLin, Ntotal]=size(Raw_PC);

indx = 1:(Nfreq-1);

ROp_data=squeeze(Raw_PC(:,:,1,:));
ROn_data=squeeze(Raw_PC(:,:,2,:));

% the simplest model of NGC employs a linear + scalar phase offset, which
% manifests as shift in the image domain.  Convert to (hybrid) image domain here:
ROp=ffts(ifft(ffts(ROp_data,1),[],1),1);
ROn=ffts(ifft(ffts(ROn_data,1),[],1),1);
phz = comp_local_pc( ROp, ROn );
phaseCoeff = reshape( phz, [ 2 size(ROp,2) size(ROp,3) ] );

%% Hello, 
% https://www.dsprelated.com/thread/4929/correlation-for-complex-number
% @jtrantow gave a very good explanation. I will try to summarize it according to my understanding. 
% 
% The correlation is a measure of similarity between signals (vectors). If we interpret signals as vectors in the N-dimensional space, the correlation becomes simply the projection of the two vectors, as @jtrantow stated. In the case, the angle between the vectors is required. Here comes the conjugate role. Let us derive it mathematically:
% 
% Let x = a exp(j theta)
% 
%     y = b exp(j phi)
% 
% Then, the correlation is defined as
% 
% x.conj(y) = a exp(j theta) . conj(b exp(j phi))
% 
%           = a exp(j theta) . (b exp(-j phi)
% 
%           = ab cos(theta - phi)
% 
% which makes sense, since we are interested in the angle between the vectors (theta - phi) not the absolute angles (theta, phi).

%% Yang's phasecorrection
% sz=size(Raw_image);
% [Nfreq, Ncoil, NphaseLine, Ntotal]=size(Raw_image);
% 
% Raw_image1 = ffts(ifft(ffts(Raw_image(:,:,:,:),1),[],1),1);
% x_index = ((1:Nfreq)-Nfreq/2)';
% 
% phase_offset = (x_index*epi_pc_lin(:)' + repmat(epi_pc_con(:)',Nfreq,1))/2;
% 
% 
% for i=1:NphaseLine
%     if mod(i,2)==1
%         phase_cor(:,:,i,:)=reshape(phase_offset,[Nfreq Ncoil Ntotal]);
%     else
%         phase_cor(:,:,i,:)=reshape(- phase_offset,[Nfreq Ncoil Ntotal]);
%     end
% end
% 
% Raw_image2 = Raw_image1.*exp(j*phase_cor);
% Raw_image3 = iffts(fft(iffts(Raw_image2,1),[],1),1);
% Raw_image_pc=reshape(Raw_image3,sz);



%% Scott's phasecorrection

sz=size(Raw_image);
[Nfreq, Ncoil, NphaseLine, Ntotal]=size(Raw_image);
k_a=zeros(NphaseLine, Nfreq, Ncoil, Ntotal);
k_b=zeros(NphaseLine, Nfreq, Ncoil, Ntotal);

k_image=permute(Raw_image(:,:,:,:),[3 1 2 4]);

k_a(1:2:end,:,:,:)=k_image(1:2:end,:,:,:);
k_b(2:2:end,:,:,:)=k_image(2:2:end,:,:,:);
[a_out,b_out] = phzapply2(k_a,k_b,phaseCoeff);
Raw_image_pc=reshape(permute(a_out+b_out,[2 3 1 4]),sz);



