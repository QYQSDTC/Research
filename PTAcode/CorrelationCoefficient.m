%%%% calculate the correlation coefficient
load('/Users/qianyiqian/desktop/matlabprograms/PTAcode/GWB/GWB.mat');
Np=getfield(simParams,'Np');
kp=getfield(simParams,'kp');
N=getfield(simParams,'N');
cthetaC=zeros(Np,Np);%%%cos(theta) between every two pulsar
r=zeros(Np,Np);%%%correlation coefficients
for i=1:1:Np
    for j=1:1:Np
cthetaC(i,j)=kp(i,:)*kp(j,:)';%%% cos(theta) between every two pulsar
    end
end
thetaC=acos(cthetaC);%%%% theta between every two pulsar

for i=1:1:Np
    for j=1:1:Np
    r(i,j) = timingResiduals_tmp(i,:)*timingResiduals_tmp(j,:)'*cthetaC(i,j)/N;
    end
end
plot(thetaC,r,'.m');