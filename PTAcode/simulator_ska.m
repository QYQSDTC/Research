% It is a simulator based on simulator2.m, we make it to simulate the data 
% for the SKA pulsars. RA, dec, distP are from..., sd is calculated based
% on the assumptions used in R. Smits ...
% Aug. 19, 2015

clear;
tic;


% master directory for simulated data
simDataDir = 'simDataSKA_locX';
searchParamsFile = 'searchParams_simDataSKA';
mkdir(simDataDir);

rng('shuffle')  % initialize the random number generator using a different seed

% ======  Useful constants  ======
pc2ly=3.261563777;  % 1 pc=3.26 ly (Julian)
dy2yr=1.0/365.25;  % 1 day=365.25 yr (Julian)  ??
kilo=1.0*10^3;  % kilo 1000

% ==== GW source parameters ====
Nsnr=3;  % number of snr of matched filter  100, 30, 8
snr_net=zeros(Nsnr,1);  % network snr
% snr_net(1)=8.0*12.7; % loc=3
% snr_net(2)=30.0*12.7;
% snr_net(3)=100.0*12.7;

% snr_net(1)=8.0*12.7/1.34/85.48*1.14;  % loc=3
% snr_net(2)=30.0*12.7/1.34/85.48*1.14;
% snr_net(3)=100.0*12.7/1.34/85.48*1.14;

% snr_net(1)=30.0*12.7/1.34/85.48;  % loc=4
% snr_net(2)=60.0*12.7/1.34/85.48;
% snr_net(3)=100.0*12.7/1.34/85.48;

% snr_net(1)=30.0*12.7/1.34/85.48*1.14/1.28;  % loc=5
% snr_net(2)=60.0*12.7/1.34/85.48*1.14/1.28;
% snr_net(3)=100.0*12.7/1.34/85.48*1.14/1.28;

snr_net(1)=30.0*12.7/1.34/85.48*1.14/1.58/0.873;  % loc=6
snr_net(2)=60.0*12.7/1.34/85.48*1.14/1.58/0.873;
snr_net(3)=100.0*12.7/1.34/85.48*1.14/1.58/0.873;

%snr_net(1)=8.0*12.7/1.34/1.54/85.48;  % loc=5
%snr_net(2)=30.0*12.7/1.34/1.54/85.48;
%snr_net(3)=100.0*12.7/1.34/1.54/85.48;

snr_tmp=zeros(Nsnr,1);
snr_tmp(1)=1.0/13.2242*snr_net(1); %10.0; %0.5; %8.0;
snr_tmp(2)=1.0/13.2242*snr_net(2); %10.0;
snr_tmp(3)=1.0/13.2242*snr_net(3);  %12.0;

Nloc=6;  % number of source sky location
alpha_tmp=zeros(Nloc,1);
delta_tmp=zeros(Nloc,1);
alpha_tmp(1)=pi/4+1.2;
delta_tmp(1)=pi/4-0.16;
alpha_tmp(2)=4.367;  % cond(A)=1.017
delta_tmp(2)=0.8796;
% location used in Taylor et al.
alpha_tmp(3)=0.9999;
delta_tmp(3)=0.5007;  % pi/2 - acos(0.48)
% loc 4, NOT Taylor et al.
alpha_tmp(4)=0.0;
delta_tmp(4)=1.2;
% loc 4, NOT Taylor et al.
alpha_tmp(5)=0.0;
delta_tmp(5)=0.0;
alpha_tmp(6)=3.5;
delta_tmp(6)=0.3;

Nomg=3;  % number of GW frequency
omega_tmp=zeros(Nomg,1);
omega_tmp(1)=2*pi/0.3925; %16.0081;  % 0.3925 yr  high freq
omega_tmp(2)=4.2743;  % 1.47 yr  low freq
omega_tmp(3)=2*pi/(1.6);  % orbital freq is 10^(-8) Hz

% optimized value
iota=0.4949;  %0.0;  %pi/4+0.6;  % inclination between orbital plane and plane of the sky
thetaN=0.5;  %pi/2;  % angle to the line of nodes
phi0=2.89;  %0.0;  % initial orbital phase

% ==== Constructing a pulsar timing array using Np pulsars ====
% read in the pulsar catalogue simulated for SKA
skamsp=load('/Users/qianyiqian/desktop/matlabprograms/PTAcode/survey_ska.mat');
[~,I]=sort(skamsp.D);
Np=1000;  % number of pulsars in the timing array

% sky location of the pulsars in the equatorial coordinate
% we need to transfer from hr angle and degree to radian
alphaP=zeros(Np,1);  % right ascension, in radian
deltaP=zeros(Np,1);  % declination, in radian
distP=zeros(Np,1);  % (parallax) distance from SSB to pulsars, from mas/pc to ly
kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,
sd=zeros(Np,1);  % standard deviation of noise for different pulsar

InList=zeros(1,2);
OutList=zeros(1,2);

for i=1:1:Np
    %     if OutList(i,1)>0
    %         alphaP(i)=OutList(i,1);  % in rad
    %     else
    %         alphaP(i)=OutList(i,1)+2*pi;  % in rad
    %     end
    
    InList(1,1)=skamsp.l(I(i));
    InList(1,2)=skamsp.b(I(i));
    
    % transfer the galactic coord. to equatorial coord.
    %[OutList,TotRot]=coco(InList,'g','j2000.0','d','r');
    [OutList,~]=coco(InList,'g','j2000.0','d','r');
    
    alphaP(i)=OutList(1,1);  % in rad
    deltaP(i)=OutList(1,2);  % in rad
    distP(i)=skamsp.D(I(i));  %0.28*kilo*pc2ly;  % in ly
    sd(i)=1.0*10^(-7);  %79.2*10^(-8); %0.148;
    
end

% % -------------------------------
% % transfer (hr,min,sec) and (degree,min,sec) to radian; mas/kpc to ly
% tmp1='J0030+0451';
% alphaP(1)=(0*15+30*15/60)*pi/180;
% deltaP(1)=(4+51/60)*pi/180;
% distP(1)=0.28*kilo*pc2ly;  % in ly
% sd(1)=79.2*10^(-8); %0.148;
% 
% tmp2='J0437-4715';
% alphaP(2)=(4*15+37*15/60)*pi/180;
% deltaP(2)=-(47+15/60)*pi/180;
% distP(2)=0.156*kilo*pc2ly;
% sd(2)=6.9*10^(-8); %0.178;
% 
% tmp3='J1640+2224';
% alphaP(3)=(16*15+40*15/60)*pi/180;
% deltaP(3)=(22+24/60)*pi/180;
% distP(3)=1.19*kilo*pc2ly;
% sd(3)=41.0*10^(-8); %0.409;
% 
% tmp4='J1713+0747';
% alphaP(4)=(17*15+13*15/60)*pi/180;
% deltaP(4)=(7+47/60)*pi/180;
% distP(4)=1.05*kilo*pc2ly;
% sd(4)=13.6*10^(-8); %0.03;
% 
% tmp5='J1744-1134';
% alphaP(5)=(17*15+44*15/60)*pi/180;
% deltaP(5)=-(11+34/60)*pi/180;
% distP(5)=0.42*kilo*pc2ly;
% sd(5)=36.6*10^(-8);
% 
% tmp6='J1857+0943';
% alphaP(6)=(18*15+57*15/60)*pi/180;
% deltaP(6)=(9+43/60)*pi/180;
% distP(6)=0.9*kilo*pc2ly;
% sd(6)=40.2*10^(-8); %0.787;
% 
% tmp7='J1909-3744';
% alphaP(7)=(19*15+9*15/60)*pi/180;
% deltaP(7)=-(37+44/60)*pi/180;
% distP(7)=1.26*kilo*pc2ly;
% sd(7)=10.0*10^(-8); %0.038;
% 
% tmp8='J1939+2134';
% alphaP(8)=(19*15+39*15/60)*pi/180;
% deltaP(8)=(21+34/60)*pi/180;
% distP(8)=5.0*kilo*pc2ly;
% sd(8)=14.1*10^(-8); %0.276;
% 
% tmp9='J2317+1439';
% alphaP(9)=(23*15+17*15/60)*pi/180;
% deltaP(9)=(14+39/60)*pi/180;
% distP(9)=1.89*kilo*pc2ly;
% sd(9)=41.2*10^(-8);
% 
% pname={tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9};c
tmp1='place_holder';
pname={tmp1};

%%
% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end

% -------------------------------
% starting epoch of the observations
start=53187;  % Modified Julian Day, 'July 1, 2004'
%finish=start+ceil(365.25*5);  % approximately, set 5 yrs
deltaT=14;  % observation cadence, in days, set biweekly
%N=389;  % 14.9 yr, 128;  % number of biweekly observations for all pulsars, fft prefer 2^n
N=130;  % 5 yrs biweekly
%N=260;  % 10 yrs biweekly
dy=zeros(N,1);  % observation epoch, in day
yr=zeros(N,1);  % observation epoch, in year
for i=1:1:N
    dy(i)=start+(i-1)*deltaT;  % dates conducting observations, in MJD
    yr(i)=2004.5+(i-1)*deltaT*dy2yr;
end

simParams = struct('Np',Np,'N',N,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp);


% -------------------------------
% calculate timing residuals induced by GW for each pulsar
%Amp=1.3*10^(-7);  % overall amplitude of timing residuals, the same for all pulsars, sec
timingResiduals=zeros(Np,N);  % signal, i.e. GW induced timing residuals, Np pulsars, N observations
timingResiduals_tmp=zeros(Np,N);   % signal without noise
phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance

Nrlz=50;  %500;  %200;  %25;  % number of noise realizations H1
noise=zeros(Np,N);  % noise
%sd=0.1*10^(-7);  % standard deviation of the normal distribtion (sec)

Nnis=50; %500; %100;  % number of realization of noise only cases H0

CA4filenames=cell(Nrlz*3+Nnis,1);  % cell array for file names of simulated data

perfect_fitness=0.0;  % fitness value for the true parameters

% set the range of the parameters
xmaxmin=zeros(7,2);  % x_max, x_min for each parameter x
xmaxmin(1,1)=2*pi;  % alpha
xmaxmin(1,2)=0.0;
xmaxmin(2,1)=pi/2;  % delta 
xmaxmin(2,2)=-pi/2;
xmaxmin(3,1)=20.0;  % angular velocity -- omega for GW
xmaxmin(3,2)=2.0;
xmaxmin(4,1)=pi;  % initial phase
xmaxmin(4,2)=0;
xmaxmin(5,1)=-8.0;  %-5.0;  %10^(-6);  % amplitude, in sec
xmaxmin(5,2)=-12.0;  % -10.0;  %10^(-8);
xmaxmin(6,1)=pi;  % inclination
xmaxmin(6,2)=0;
xmaxmin(7,1)=pi;  % polarization
xmaxmin(7,2)=0;

save(searchParamsFile,'xmaxmin');

%Standardized true parameter values
%stdTrueCoord = zeros(1,12);  % 4+Np 
stdTrueCoord = zeros(1,7); % 7 parameters other than pulsar phases

% calculate SNR
%snr_chr=0;  % characteristic srn = <\rho> = sqrt[(h,h)] srn of signal -- the strenghth of signal
%snr_mf=0;  % snr of matched filter = \rho = x/\sigma, normalized matched filter

nf=0;  % counter for number of files
for ii=1:1:Nsnr
    
    genHypothesis='H1 data';
    Amp=snr_tmp(ii)*2*10^(-9);
    
    for j=6:1:6  %Nloc
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta_tmp(j))*cos(alpha_tmp(j));
        k(2)=cos(delta_tmp(j))*sin(alpha_tmp(j));
        k(3)=sin(delta_tmp(j));
        
        for l=3:1:3  %Nomg
            
            snr_id=ii; %num2str(ii);
            loc_id=j; %num2str(j);
            omg_id=l; %num2str(l);
            
            snr_chr=0.0;  % initialize to zero
            
            % simulate the signals
            for i=1:1:Np
                
                theta=acos(k*kp(i,:)');
                %sprintf('%d pulsar theta=%g',i,theta)
                %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
                %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
                phiI(i)=mod(phi0-0.5*omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
                
                %disp(['pulsar = ', num2str(i), ' ', num2str(phiI(i))])
                
                timingResiduals_tmp(i,:)=FullResiduals(alpha_tmp(j),delta_tmp(j),omega_tmp(l),phi0,phiI(i),alphaP(i),deltaP(i),...
                    Amp,iota,thetaN,theta,yr);
                
                %fftsignal(i,:)=fft(timingResiduals_tmp(i,:));

                % calculate the perfect fitness value
                
                snr_chr = snr_chr + dot(timingResiduals_tmp(i,:),timingResiduals_tmp(i,:)) / sd(i)^2;
                
                % standardization of the true coordinates
                stdTrueCoord(1)=(alpha_tmp(j)-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
                stdTrueCoord(2)=(delta_tmp(j)-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
                stdTrueCoord(3)=(omega_tmp(l)-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));  % [2, 20]
                stdTrueCoord(4)= mod(phi0,pi)/pi;  % [0, pi]
                stdTrueCoord(5)=(log10(Amp)-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
                stdTrueCoord(6)=(iota-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
                stdTrueCoord(7)=(thetaN-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));
                
            end
            
            
            %snr_chr=sqrt(snr_chr/Np);  % averaged snr--root of mean square of individual snr
            snr_chr=sqrt(snr_chr);
            
            % signal + noise realizations
            for jj=1:1:Nrlz
                
                nf=nf+1;
                
                rlz_id=jj; %num2str(jj);
                
                % structure to store the id tag for each metadata file
                id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
                
                %snr_tmp=0;
                
                for i=1:1:Np
                    
                    % generating a realization of noise
                    noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
                    % calculate the actual snr
                    %fftnoise(i,:)=fft(noise(i,:));
                    
                    %snr_tmp=snr_tmp+ fftsignal(i,:)*fftsignal(i,:)/fftnoise(i,:);
                    timingResiduals(i,:)=timingResiduals_tmp(i,:)+noise(i,:);  % add noise on signal
                    
                end
                
                inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
                    'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,...
                    'xmaxmin',xmaxmin);
                
                perfect_fitness = LLR_PSOmpp(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
                %true_fitness = fitnessTrue_ie(alpha_tmp(j),delta_tmp(j),omega_tmp(l),Amp,iota,thetaN,phi0,phiI,inParams);
                
                %disp(['In simulator2: perfect_fitness: ', num2str(perfect_fitness)]);
                
                % save metadata into a file for each realization (file name rule)
                filename=strcat('snr',num2str(ii),'loc',num2str(j),'omg',num2str(l),'rlz',num2str(jj),'.mat');
                snr=snr_tmp(ii);
                alpha=alpha_tmp(j);
                delta=delta_tmp(j);
                omega=omega_tmp(l);
                CA4filenames{nf}=filename;
                
                save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
                    'timingResiduals','noise','yr', 'pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');
                             
            end
            
        end
        
    end
    
end


% simulating noise without signal
%Nnis=10; %100;  % number of realization of noise
perfect_fitness = 0;
snr_chr = 0;
timingResiduals_tmp=0;
for ii=1:1:Nnis
    
    genHypothesis='H0 data';
    snr_id=0; %'0';
    loc_id=0; %'0';
    omg_id=0; %'0';
    rlz_id=ii; %num2str(ii);
    
    id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);

    for i=1:1:Np
        
        % generating a realization of noise
        noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
        timingResiduals(i,:)=noise(i,:);
        %timingResiduals(i,:)=timingResiduals_tmp(i,:)+noise(i,:);  % add noise on signal
        
    end
    
    
    % save metadata into a file for each realization (file name rule)
    filename=strcat('noise',num2str(ii),'.mat');
    snr=0.0;  %snr_tmp(ii);
    alpha=0.0;  %alpha_tmp(j);
    delta=0.0;  %delta_tmp(j);
    omega=0.0;  %omega_tmp(l);
    iota=0.0;
    thetaN=0.0;
    phi0=0.0;
    save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
        'timingResiduals','noise','yr','pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp');
    CA4filenames{nf+ii}=filename;
    
end


toc; % stop the stopwatch timer
% end of script
