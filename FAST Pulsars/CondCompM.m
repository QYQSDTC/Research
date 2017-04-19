% calculate the condition number for the companion matrix D used for
% solving the quartic equaitons in the maximization of pulsar phases
% algorithm
% Oct. 29, 2014.  need to change Dir and Amp

clear;

tic;

filename='snr3loc3omg3rlz10.mat';
% load the data into working space
simDataDir='/home/juan/research/PULSARTIMING/MultiCW/simData9_snr123_loc3_sr/';
resDataDir='/home/juan/research/PULSARTIMING/MultiCW/simData9_snr123_loc3_sr/results/';

data = load([simDataDir,filesep,filename]);
load([resDataDir,filesep,filename],'bestLocation','fitnessVals');

% set the range of the parameters
tmpmaxmin = load('/home/juan/research/PULSARTIMING/MultiCW/searchParams_simDataX.mat');
xmaxmin = tmpmaxmin.xmaxmin;


Np=data.simParams.Np;  % number of pulsar in PTA
N=data.simParams.N;
alphaP=data.simParams.alphaP;
deltaP=data.simParams.deltaP;
%MulResiduals=PTAParams.MulResiduals;  % signals from individual sources
timingResiduals_noise=data.timingResiduals;  % signals with noise
kp=data.simParams.kp;
sd=data.simParams.sd;
yr=data.yr;
%phiI=data.simParams.phiI;

% PSO input parameters
inParams = struct('Np',Np,'N',N,'s',timingResiduals_noise,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,...
                  'xmaxmin',xmaxmin);

[fitVal_bf,allParam] =  LLR_PSOmpp(bestLocation,inParams);  % - LogLikelihoodRatio, minimization
%bfResiduals=bfResidualsfunc(allParam,inParams);

% parameter values from bestLocation
omega=allParam(3);
phi0=allParam(4);
Amp=allParam(5);
iota=allParam(6);
thetaN=allParam(7);
%phiI=zeros(Np,1);
%phiI=allParam(8:length(allParam));

% parameter values from true signal
omega_t=data.omega;
phi0_t=data.phi0;
%Amp_t=allParam(5);
iota_t=data.iota;
thetaN_t=data.thetaN;
% choose correct Amp
%Amp_t=1.5124e-08/1.8;   % snr3
%Amp_t=1.5124e-09/1.8;   % snr2
%Amp_t=7.5619e-10;   % snr1

% simData9
Amp_t=1.9207e-07;  % snr3 net_snr=100
%Amp_t=1.9207e-07*3.0/10.0;  % snr2 net_snr=30
%Amp_t=1.9207e-07*0.8/10.0;  % snr1 net_snr=8

% searching ra and dec
Na=50; %100;
Nd=25; %50;
alpha=zeros(Na,1);  % ra of source
delta=zeros(Nd,1);  % dec of source
%ks=zeros(1,3);
D=zeros(4,4);
Dcond=zeros(Na,Nd,Np);  % condition number of D at a sky location
D_t=zeros(4,4);
Dcond_t=zeros(Na,Nd,Np);  % condition number of D at a sky location
%Inten=zeros(Na,Nd);

stdTrueCoord = zeros(1,7);  % std coord. for LogLikelihoodRatioMP.m
stdTrueCoord(3)=(omega-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));
stdTrueCoord(4)= mod(phi0,pi)/pi;  % [0, pi]
stdTrueCoord(5)=(log10(Amp)-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
stdTrueCoord(6)=(iota-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
stdTrueCoord(7)=(thetaN-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));

stdTrueCoord_t = zeros(1,7);  % std coord. for LogLikelihoodRatioMP.m
stdTrueCoord_t(3)=(omega_t-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));
stdTrueCoord_t(4)= mod(phi0_t,pi)/pi;  % [0, pi]
stdTrueCoord_t(5)=(log10(Amp_t)-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
stdTrueCoord_t(6)=(iota_t-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
stdTrueCoord_t(7)=(thetaN_t-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));

% coefficients for the characteristic polynomial
e=zeros(Np,5);
e_t=zeros(Np,5);

% calculate the condition number for the companion matrix D
for i=1:1:Na
    alpha(i)=(i-1)*2*pi/Na;
    
    for j=1:1:Nd
        delta(j)=(j-1)*pi/Nd-pi/2;
        
        % calculate the fitness at true intrinsic params
        stdTrueCoord(1)= (alpha(i)-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
        stdTrueCoord(2)= (delta(j)-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
 
        stdTrueCoord_t(1)=stdTrueCoord(1);
        stdTrueCoord_t(2)=stdTrueCoord(2);
        
        for k=1:1:Np
            
            % for bestLocation
            D=zeros(4,4);  % reset D to be zero matrix for each pulsar
            D(2,1)=1.0; % initialize D
            D(3,2)=1.0;
            D(4,3)=1.0;
            
            e = quarticCoeff(stdTrueCoord,inParams);
            D(1,1)=-e(k,2)/e(k,1);
            D(1,2)=-e(k,3)/e(k,1);
            D(1,3)=-e(k,4)/e(k,1);
            D(1,4)=-e(k,5)/e(k,1);
            
            Dcond(i,j,k)=cond(D,2);  % cond(D,inf)
            
            % for true location
            D_t=zeros(4,4);  % reset D to be zero matrix for each pulsar
            D_t(2,1)=1.0; % initialize D
            D_t(3,2)=1.0;
            D_t(4,3)=1.0;
            
            e_t = quarticCoeff(stdTrueCoord_t,inParams);
            D_t(1,1)=-e_t(k,2)/e_t(k,1);
            D_t(1,2)=-e_t(k,3)/e_t(k,1);
            D_t(1,3)=-e_t(k,4)/e_t(k,1);
            D_t(1,4)=-e_t(k,5)/e_t(k,1);
            
            Dcond_t(i,j,k)=cond(D_t,2);  % cond(D,inf)
                           
        end    
        
    end
end

% find the maximum of the conditional number among pulsars
Dcond_max=zeros(Na,Nd);
Dcond_t_max=zeros(Na,Nd);

for i=1:1:Na
    for j=1:1:Nd
        
        Dcond_max(i,j)=max(Dcond(i,j,:));
        Dcond_t_max(i,j)=max(Dcond_t(i,j,:));
        
    end
end


% plot the condition numbers for D 
figure  % bestLocation
surf(alpha,delta,log10(Dcond_max)');
xlim([0 2*pi]);
ylim([-pi/2 pi/2]);
xlabel('Right Ascension \alpha');
ylabel('Declination \delta');
title('Condition number of matrix D at different sky location and pso parameters');

[v,ind]=max(log10(Dcond_max));
[v1,ind1]=max(max(log10(Dcond_max)));
disp(['The largest element in Dcond_max is: ' num2str(v1) ' at (' num2str(ind(ind1)) ',' num2str(ind1) ')']);
disp(['The values for ra and dec are: ', num2str(alpha(ind(ind1))),' and ',num2str(delta(ind1))]);


figure  % true
surf(alpha,delta,log10(Dcond_t_max)');
xlim([0 2*pi]);
ylim([-pi/2 pi/2]);
xlabel('Right Ascension \alpha');
ylabel('Declination \delta');
title('Condition number of matrix D at different sky location and true parameters');

[w,ind]=max(log10(Dcond_t_max));
[w1,ind1]=max(max(log10(Dcond_t_max)));
disp(['The largest element in Dcond_t_max is: ' num2str(w1) ' at (' num2str(ind(ind1)) ',' num2str(ind1) ')']);
disp(['The values for ra and dec are: ', num2str(alpha(ind(ind1))),' and ',num2str(delta(ind1))]);



toc;

% End