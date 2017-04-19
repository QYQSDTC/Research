%% stochastic gravitational background
function [timingResidulas_tmp]=SGWB(alphaP,deltaP,distP,yr,Ns,kp)
        
[Amp,alpha_tmp,delta_tmp,fgw,iota,thetaN,phi0,r]=GenerateRandomGWSource(Ns);

omega_tmp = 2 .* pi .* fgw .* 3 .* 10^7;%% change the unit to yr^-1
timingResiduals_tmp=zeros(Np,N);   % signal without noise
phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance
tmp=zeros(1,N);
for i=1:1:Np
    for j=1:1:Ns  % number of GW sources
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta_tmp(j))*cos(alpha_tmp(j));
        k(2)=cos(delta_tmp(j))*sin(alpha_tmp(j));
        k(3)=sin(delta_tmp(j));
        theta=acos(k*kp(i,:)');
        %sprintf('%d pulsar theta=%g',i,theta)
        %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
        %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
        phiI(i)=mod(phi0(j)-0.5*omega_tmp(j)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        %disp(['pulsar = ', num2str(i), ' ', num2str(phiI(i))])
        
        tmp = FullResiduals(alpha_tmp(j),delta_tmp(j),omega_tmp(j),phi0(j),phiI(i),alphaP(i),deltaP(i),...
            Amp(j),iota(j),thetaN(j),theta,yr);
        timingResiduals_tmp(i,:) = timingResiduals_tmp(i,:)+tmp';
        %fftsignal(i,:)=fft(timingResiduals_tmp(i,:));
        
        % calculate the perfect fitness value
        
        %snr_chr = snr_chr + dot(timingResiduals_tmp(i,:),timingResiduals_tmp(i,:)) / sd(i)^2;
        
        % standardization of the true coordinates
        %                 stdTrueCoord(1)=(alpha_tmp(j)-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
        %                 stdTrueCoord(2)=(delta_tmp(j)-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
        %                 stdTrueCoord(3)=(omega_tmp(j)-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));  % [2, 20]
        %                 stdTrueCoord(4)= mod(phi0(j),pi)/pi;  % [0, pi]
        %                 stdTrueCoord(5)=(log10(Amp(j))-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
        %                 stdTrueCoord(6)=(iota(j)-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
        %                 stdTrueCoord(7)=(thetaN-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));
        
    end

    
    %snr_chr=sqrt(snr_chr/Np);  % averaged snr--root of mean square of individual snr
    %snr_chr=sqrt(snr_chr);
    
    % signal + noise realizations
    %for jj=1:1:Nrlz
    
    %nf=nf+1;
    
    %rlz_id=jj; %num2str(jj);
    
    % structure to store the id tag for each metadata file
    %id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
    
    %snr_tmp=0;
    
    %                 for i=1:1:Np
    %
    %                     % generating a realization of noise
    %                     noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
    %                     % calculate the actual snr
    %                     %fftnoise(i,:)=fft(noise(i,:));
    %
    %                     %snr_tmp=snr_tmp+ fftsignal(i,:)*fftsignal(i,:)/fftnoise(i,:);
    %                     timingResiduals(i,:)=timingResiduals_tmp(i,:)+noise(i,:);  % add noise on signal
    %
    %                 end
    
%     inParams = struct('Np',Np,'N',N,'Ns',Ns,'s',timingResiduals_tmp,'sd',sd,...
%         'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr);
    
    %perfect_fitness = LLR_PSOmpp(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
    %true_fitness = fitnessTrue_ie(alpha_tmp(j),delta_tmp(j),omega_tmp(l),Amp,iota,thetaN,phi0,phiI,inParams);
    
    %disp(['In simulator2: perfect_fitness: ', num2str(perfect_fitness)]);
    
    % save metadata into a file for each realization (file name rule)
    %filename=strcat('snr',num2str(ii),'loc',num2str(j),'omg',num2str(j),'rlz',num2str(jj),'.mat');
    
    %end
    
    % end
    
end