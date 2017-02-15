%% Main function of the algorithm that maximizes pulsar phases
% Based on the note for the maximization of pulsar phases in a pulsar timing array
% Yan Wang and Soumya Mohanty, UT-Brownsville, May 26, 2014
% One GW source case. PSO script is added.

tic;

% clear all;
%clc;

%% function generating the simulated timing residuals induced by multiple 
% GW sources for each pulsar in a PTA, source parameters and pulsar
% parameters should be set in the function 'MulPTAsimulator'
[SrcParams, PTAParams, pname] = MulPTAsimulator;

%% transfer the data from 'MulPTAsimulator' into Main
Ns=SrcParams.Ns;  % number of GW source
alpha=SrcParams.alpha;
delta=SrcParams.delta;
omega=SrcParams.omega;
phi0=SrcParams.phi0;
Amp=SrcParams.Amp;
iota=SrcParams.iota;
thetaN=SrcParams.thetaN;

Np=PTAParams.Np;  % number of pulsar in PTA
N=PTAParams.N;
alphaP=PTAParams.alphaP;
deltaP=PTAParams.deltaP;
MulResiduals=PTAParams.MulResiduals;  % signals from individual sources
MulResiduals_noise=PTAParams.MulResiduals_noise;  % signals with noise
kp=PTAParams.kp;
sd=PTAParams.sd;
yr=PTAParams.yr;
phiI=PTAParams.phiI;
%pname=PTAParams.pname;  % transfer via output parameter of a function

disp('True Intrinsic Parameters in Real Coord.')
disp([alpha,delta,omega,phi0,Amp,iota,thetaN]');
disp('True Extrinsic Parameters in Real Coord.')
disp(phiI);

%% set the range of the parameters
xmaxmin=zeros(7,2);  % x_max, x_min for each parameter x
xmaxmin(1,1)=2*pi;  % alpha
xmaxmin(1,2)=0;
xmaxmin(2,1)=pi/2;  % delta 
xmaxmin(2,2)=-pi/2;
xmaxmin(3,1)=20.0;  % angular velocity -- omega for GW
xmaxmin(3,2)=2.0;
xmaxmin(4,1)=pi;  % initial phase
xmaxmin(4,2)=0.0;
xmaxmin(5,1)=-6.0;  %10^(-6);  % amplitude, in sec
xmaxmin(5,2)=-8.0;  %10^(-8);
xmaxmin(6,1)=pi;  % inclination
xmaxmin(6,2)=0;
xmaxmin(7,1)=pi;  % polarization
xmaxmin(7,2)=0;

% Standardized [0,1] true parameter values for intrinsic parameters
stdTrueCoord = zeros(1,7); 
stdTrueCoord(1)= (alpha-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
stdTrueCoord(2)=(delta-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
stdTrueCoord(3)=(omega-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));
stdTrueCoord(4)= mod(phi0,pi)/pi;  % [0, pi]
stdTrueCoord(5)=(log10(Amp)-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
stdTrueCoord(6)=(iota-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
stdTrueCoord(7)=(thetaN-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));

%% timing residuals are the sum of contributions from all individual sources
timingResiduals=zeros(Np,N);  % noise free timing residuals
timingResiduals_noise=zeros(Np,N);  % timing residuals contain noise
for j=1:1:Np
    for i=1:1:Ns
        
        timingResiduals_noise(j,:)=timingResiduals_noise(j,:)+MulResiduals_noise(j,:,i);
        % CAUSION: noise should be added once for each pulsar (note change code in MulPTAsimulator)
        timingResiduals(j,:)=timingResiduals(j,:)+MulResiduals(j,:,i);
        
    end
end

% PSO input parameters
inParams = struct('Np',Np,'N',N,'s',timingResiduals_noise,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,...
                  'xmaxmin',xmaxmin);

% plot the timing residuals from 'MulPTAsimulator'
figure
for i=1:1:Np
    subplot(4,2,i)
    %plot(yr,MulResiduals(i,:),'.-');
    plot(yr,timingResiduals(i,:),'.-');
    grid on;
    hold on
    plot(yr,timingResiduals_noise(i,:),'r.-');
    xlabel('years');
    ylabel('residuals (sec)');
    title(pname(i));
end

%% test the functional form for the maximum of phase method
timingResiduals_test=zeros(Np,N);  % from formular 1
timingResiduals_test2=zeros(Np,N);  % from formular 2

% sky location in Cartesian coordinate
k=zeros(Ns,3);  % unit vector pointing from SSB to source
for i=1:1:Ns
    k(i,1)=cos(delta(i))*cos(alpha(i));
    k(i,2)=cos(delta(i))*sin(alpha(i));
    k(i,3)=sin(delta(i));
end

Phi=omega*yr;  % Phi=omega*t, N by 1 matrix
for i=1:1:Np
   
    theta=acos( dot(k(1,:),kp(i,:)) );  % only one source here
    alphatilde=alpha(1)-alphaP(i);
    
    Pp=-cos(deltaP(i))^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta(1))^2)...
        +sin(deltaP(i))^2*cos(delta(1))^2-0.5*sin(2*deltaP(i))*cos(alphatilde)*sin(2*delta(1));
    
    Pc=2*cos(deltaP(i))*sin(alphatilde)*(cos(deltaP(i))*cos(alphatilde)*sin(delta(1))...
        -sin(deltaP(i))*cos(delta(1)));
    
    Fp=Pp/(1-cos(theta));
    Fc=Pc/(1-cos(theta));
    
    %disp(['in testMul','Pp',num2str(Pp),'Pc',num2str(Pc),'theta',num2str(theta)])
    
    A=2*Amp*sqrt( (1+cos(iota)^2)^2*(Fp*cos(2*thetaN)-Fc*sin(2*thetaN))^2 ...
        + 4*cos(iota)^2*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))^2 ); % A is different for pulsars
    
    tmp=-2*cos(iota)/(1+cos(iota)^2)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))/(Fp*cos(2*thetaN)-Fc*sin(2*thetaN));
    % solve psi, atan or atan2 ?
    %psi=atan(tmp);
    psi=atan2( -2*cos(iota)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN)),(1+cos(iota)^2)*(Fp*cos(2*thetaN)-Fc*sin(2*thetaN)) );
    
    sin_ts=sin(phi0+psi+Phi);  % N by 1 matrix
    cos_ts=cos(phi0+psi+Phi);
    
    B=0.5*A*sin(phi0)*sin_ts;
    C=-0.5*A*cos(phi0)*sin_ts;
    D=0.5*A*sin(phi0)*cos_ts;
    E=-0.5*A*cos(phi0)*cos_ts;

    timingResiduals_test(i,:)=(B-E)*cos(2*phiI(i))+(C+D)*sin(2*phiI(i))+(B+E);  % Eq 12
    timingResiduals_test2(i,:)=A*sin(phi0-phiI(i))*sin(phi0+phiI(i)+psi+Phi);  % Eq 9

end

% figure
% for i=1:1:Np
%     subplot(4,2,i)
%     plot(yr,timingResiduals_test(i,:),'.-');
%     grid on;
% %     hold on
% %     plot(yr,bestfitResiduals(i,:),'r.-');
%     xlabel('years');
%     ylabel('residuals (sec)');
%     title(pname(i));
% end
% 
% figure
% for i=1:1:Np
%     subplot(4,2,i)
%     plot(yr,timingResiduals_test2(i,:),'.-');
%     grid on;
% %     hold on
% %     plot(yr,bestfitResiduals(i,:),'r.-');
%     xlabel('years');
%     ylabel('residuals (sec)');
%     title(pname(i));
% end


%% LLR_PSO2 accepts only the standard coordinates [0,1] as the input variable
% Here, the intrinsic parameters are: alpha, delta, omega, phi0, Amp, iota, thetaN
[fitVal, allParam_in] = LLR_PSO2(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization

% fitness for the true parameters (both intrinsic and extrinsic)
fitVal_ie = fitnessTrue_ie(alpha,delta,omega,Amp,iota,thetaN,phi0,phiI,inParams);

disp(['Fitness for true int & ext parameters: ', num2str(fitVal_ie)]);  % fitVal with true intrinsic parameters
disp(['Fitness for true intrinsic parameters: ', num2str(fitVal)]);  % fitVal with true intrinsic parameters
disp('Reconstructed ext parameters in true Coord: ')
disp(allParam_in(8:length(allParam_in)));

%% GLOBAL OPTIMIZATION using Matlab code
% fHandle = @(x) LLR_PSO2(x,inParams);
% 
% % Create a solver object
% ms = MultiStart('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxTime',3600,...
%     'UseParallel','always');
% %gs = GlobalSearch('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxTime',1000);
% 
% % Set Start Points for MultiStart (MultiStart only)
% xstart=rand(1,7);
% 
% % Create a problem structure
% opts = optimset('Algorithm','interior-point'); % Create an options structure 
% problem = createOptimProblem('fmincon', ... 
%     'x0', xstart, ...
%     'objective',fHandle,...
%     'lb',[ 0., 0., 0., 0., 0., 0., 0. ], ...
%     'ub',[ 1., 1., 1., 1., 1., 1., 1. ], ...
%     'options',opts);
% 
% % Run the solver
% instance=1000;  % number of runs
% disp('running ms ... ');
% [xmin,fmin,flag,outpt,manymins] = run(ms,problem,instance);
% %[xmin,fmin,flag,outpt,manymins] = run(gs,problem);


%% PARTICLE SWARM OPTIMIZATION, Soumya's scripts
nRuns = 1;
fHandle = @(x) LLR_PSO2(x,inParams);
nDim = length(xmaxmin(:,1));
% PSO configuration parameter structure
P = psoparamstruct(1,'default');
P.tuningVars.numPart = 40;
P.tuningVars.maxSteps = 800;
P.tuningVars.inertiaDecayLaw = 'linear';
P.tuningVars.inertiaDecayParams = [0.9,0.4,P.tuningVars.maxSteps,0.2];
P.topology.topoScheme = 'ring';
P.topology.topoParams = 3;
P.convergeScheme = 'maxSteps';
% PSO output parameter structure
outP = struct('outFileNames',{{'log',''}},'graphics','off','status','off');
 %Start independent PSO runs
 bestLocationVec = zeros(nRuns,nDim);
 %bestLocRealC = zeros(nRuns,nDim);
 bestFitValVec = zeros(nRuns,1);
 nFuncEvalsVec = zeros(nRuns,1);
 nIterVec = zeros(nRuns,1);
 wallClkTimeVec = zeros(nRuns,1);
    %Log info
%     fprintf(fidLog,'\t fitness found/number of function evaluations \n ----\n');
%     disp('fitness found/number of function evaluations');
for lprun = 1:nRuns
    
    tic;
    psoResults=pso(fHandle,nDim,P,outP);
    wallClkTime = toc;
    wallClkTimeVec(lprun) = wallClkTime/60;%min
    
    bestLocationVec(lprun,:)=psoResults.bestLocation;
    %FIXME the line below has a bug
    %bestLocRealC(lprun,:) = LLR_PSO2(psoResults.bestLocation,inParams);
    bestFitValVec(lprun) = psoResults.bestSNR;
    disp(['Fitness value from run ',num2str(lprun),': ',num2str(bestFitValVec(lprun))]);
    nFuncEvalsVec(lprun) = psoResults.totalFuncEvals;
    nIterVec(lprun) = psoResults.totalSteps;
    
    %Log info
    %         fprintf(fidLog,'\t %e/', bestFitValVec(lprun));
    %         fprintf(fidLog,'%u', ...
    %             nFuncEvalsVec(lprun));
end

[bestFitVal,bestFitIndx] = min(bestFitValVec);
%     fprintf(fidLog,'\n\t Best Run: %d \n',bestFitIndx);
bestLocation = bestLocationVec(bestFitIndx,:);

%[~,realC]=LLR_PSO2(bestLocation,inParams);

%% calculate and plot the best-fit residuals comparing with the true residuals
[fitVal_bf,allParam] = LLR_PSO2(bestLocation,inParams);
%[fitVal,allParam] = LLR_PSO2(intPstd,inParams);  % allParam intrinsic and extrinsic
%allParam(8:length(allParam))=phiI;
bfResiduals=bfResidualsfunc(allParam,inParams);
disp(['Fitness value from LLR_PSO2 ',num2str(fitVal_bf)]);
disp(['Estimated int parameters ',num2str(allParam(1:7))]);
figure
for i=1:1:Np
    subplot(4,2,i)
    plot(yr,timingResiduals(i,:),'.-');  % true residuals in blue
    grid on;
    hold on
    plot(yr,bfResiduals(i,:),'r.-');  % best-fit residuals in red
    xlabel('years');
    ylabel('residuals (sec)');
    title(pname(i));
end


toc;
% END of MAIN