%% FAST ms Pulsars coordinates transformation
%% 2017.4.13 QYQ

fastmsp1=load('/Users/qianyiqian/desktop/matlabprograms/FAST Pulsars/fastmsp1.mat');
n = fastmsp1.X;
Np = length(n);
InList = zeros(1,2);
OutList = zeros(1,2);
alphaP = zeros(Np,1);
deltaP = zeros(Np,1);
for i = 1:1:Np
    
    InList(1,1) = fastmsp1.l(i);
    InList(1,2) = fastmsp1.b(i);
    [OutList,~]=coco(InList,'g','j2000.0','d','r');
    alphaP(i)=OutList(1,1);
    deltaP(i)=OutList(1,2);
end