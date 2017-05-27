%% Check the relationship between pulsars and condition number
%% 2017.5.27 QYQ
%%
clear;

fastmsp=load('/Users/qianyiqian/desktop/matlabprograms/FAST Pulsars/fastmsp2.mat');
[~,I]=sort(fastmsp.D);
for Np = 50:50:450  %check how every add 50 pulsars to affect the condition number 
%Np = length(fastmsp.D);
InList = zeros(1,2);
OutList = zeros(1,2);
alphaP = zeros(Np,1);
deltaP = zeros(Np,1);
kp = zeros(Np,3);
for i = 1:1:Np
    
    InList(1,1) = fastmsp.l(I(i));
    InList(1,2) = fastmsp.b(I(i));
    [OutList,~]=coco(InList,'g','j2000.0','d','r');
    alphaP(i)=OutList(1,1);
    if alphaP(i) <= 0
      alphaP(i) = 2*pi + alphaP(i);
    end
    deltaP(i)=OutList(1,2);
end
% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end
filename = strcat('FASTband2-',num2str(Np),'.mat');
save(filename,'alphaP','deltaP','kp');
end