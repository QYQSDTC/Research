%%%% calculate the correlation coefficient
function [CE,thetaC] = CorrelationCoefficient()
run simulator4.m;
load('/Users/qianyiqian/desktop/matlabprograms/PTAcode/GWB/GWB.mat');
Np=getfield(simParams,'Np');
kp=getfield(simParams,'kp');
%N=getfield(simParams,'N');
cthetaC=zeros(1,(Np-1)*Np/2);%%%cos(theta) between every two pulsar
CE=zeros(1,(Np-1)*Np/2);%%%correlation coefficients
ct = 1;%% counter

for i=1:1:Np-1
    for j=1+i:1:Np
        
        cthetaC(:,ct)=kp(i,:)*kp(j,:)';%%% cos(theta) between every two pulsar
        R=(timingResiduals_tmp(i,:)*timingResiduals_tmp(i,:)')...
       *(timingResiduals_tmp(j,:)*timingResiduals_tmp(j,:)');
        CE(:,ct) = timingResiduals_tmp(i,:)*timingResiduals_tmp(j,:)'/sqrt(R);
        ct = ct+1;
        
    end
end
thetaC=acos(cthetaC)*180/pi;%%%% theta between every two pulsar
%plot(thetaC,CE,'.k');