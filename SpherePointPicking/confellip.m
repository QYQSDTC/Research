% Calculate the confidence elliptic for a 2-D scatter plot, use KDE
% function 'kde2d' and determine the contour level to use
% Sep 26, 2016. YW
% Jan 18, 2017. Yi-Qian Qian added contour and Monte Carlo simulation 

function [SD1,SD2]=confellip(XX,YY,varargin)
if nargin < 3
   % fprintf('Use the default confidence levels 68%% and 95%%');
    disp('Use the default confidence levels 68% and 95%');
end
if nargin == 3 || nargin == 5
    error('The input argument is wrong');
end
conf1 = 0.68; %default confidence level 1
conf2 = 0.95; %default confidence level 2
pp = varargin;
while length(pp) >= 2
    prop = pp{1};
    if isa(pp{2},'double')
        val = pp{2};
    else
        error('The input argument is wrong');
    end
    pp = pp(3:end);
    switch prop
        case 'conf1'
            conf1 = val;
        case 'conf2'
            conf2 = val;
        otherwise
            error('The input argument is wrong');
    end
end
    
% load data into workspace
%results=load('/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc6/resultsSKA_snr0123_AnglesWrapped_v1_May26_2016/summary/H1stats_snr1_loc6_omg3.mat');
%results=load('F:\GitProject\MatlabPrograms\SpherePointPicking/H1stats_snr1_loc9_omg3.mat');
%results=load('/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc9/resultsSKA_snr0123_loc9/summary/H1stats_snr1_loc9_omg3.mat');

%data=[results.estSigParams(:,1), results.estSigParams(:,2)]; % x, y coords 
data(:,1)=XX;
data(:,2)=YY;
Ns=length(XX);  % sample size
Conf1 = floor(Ns*conf1);  % contour level
Conf2 = floor(Ns*conf2);  % contour level
ra=mean(XX);  % mean value of RA
dec=mean(YY);  % mean value of Dec

% calculate the density by 'kde2d' function
[bandwidth,density,X,Y]=kde2d(data);

% plot the contour to have the coordinates of different contour levels
Ncl=100;  % number of contour level
figure;
[C,h]=contour3(X,Y,density,Ncl);
level=zeros(Ncl,1);  % contour level

Np=zeros(Ncl,1);  % number of date points located inside that contour
NpConf1=zeros(3,1);  % number of random sphere points located inside that contour
NpConf2=zeros(3,1);  % number of random sphere points located inside that contour
cont=1;  % counter
xq=data(:,1);  % x coord.
yq=data(:,2);  % y coord.

% generate n random uniform points on the sphere shell use
% 'SpherePointPicking', RA=[0 2*pi], Dec=[-pi/2 pi/2]
n = 50000;  % larger number means higher precision and longer computation time
[theta,phi]=SpherePointPicking(n); 

for i=1:1:Ncl
    
    level(i)=C(1,cont);
    Nc=C(2,cont);  % number of coord. for that contour
    int=cont+1;  % initial index for the level
    fin=int+Nc-1;  % final index for the level
    cont=cont+C(2,cont)+1;
    
    %xv=zeros(Nc,1);  % x coord. of contour
    %yv=zeros(Nc,1);
    xv=C(1,int:fin)';
    yv=C(2,int:fin)';
    
    in=inpolygon(xq,yq,xv,yv);
    Np(i)=sum(in);
    
end

% generate the Conf1 contour
for i=1:1:Ncl
    if Np(i)<Conf1
        disp(['Index for ',num2str(conf1*100),'% level is: ', num2str(i), ', Np= ', num2str(Np(i))...
            ', contour level is: ', num2str(level(i))]);
        break;
    end
end

i1 = i;
cont = 1;  % reinitialize counter
for i=1:1:i1
    
    level(i)=C(1,cont);
    Nc=C(2,cont);  % number of coord. for that contour
    int=cont+1;  % initial index for the level
    fin=int+Nc-1;  % final index for the level
    cont=cont+C(2,cont)+1;
    
    %xv=zeros(Nc,1);  % x coord. of contour
    %yv=zeros(Nc,1);
    xv1=C(1,int:fin)';
    yv1=C(2,int:fin)';
    
    in=inpolygon(theta,phi,xv1,yv1);
    NpConf1(i)=sum(in);
    P = NpConf1(i)/n ; %calculate the fraction of N68 to the solid angle
end
SD1 = P*41253.0;
disp(['The fraction of ',num2str(conf1*100),'% to the solid angle is :',num2str(P)]);
disp(['Square degree for ',num2str(conf1*100),'% is :',num2str(SD1)]);
% figure;
% plot(data(:,1),data(:,2),'g+');
% hold on;
% % plot(theta,phi,'.r');
% % hold on;
% plot(xv,yv,'k');
% xlim([0 2*pi])
% ylim([-pi/2 pi/2])
% title('68% confidence contour');

% generate the N95 contour
for i=1:1:Ncl
    if Np(i)<Conf2
        disp(['Index for ',num2str(conf2*100),'% level is: ', num2str(i), ', Np= ', num2str(Np(i))...
            ', contour level is: ', num2str(level(i))]);
        break;
    end
end

i2 = i;
cont = 1;  % reinitialize counter
for i=1:1:i2
    
    level(i)=C(1,cont);
    Nc=C(2,cont);  % number of coord. for that contour
    int=cont+1;  % initial index for the level
    fin=int+Nc-1;  % final index for the level
    cont=cont+C(2,cont)+1;
    
    %xv=zeros(Nc,1);  % x coord. of contour
    %yv=zeros(Nc,1);
    xv2=C(1,int:fin)';
    yv2=C(2,int:fin)';
    
    in=inpolygon(theta,phi,xv2,yv2);
    NpConf2(i)=sum(in);
    P = NpConf2(i)/n ; % calculate the fraction of N95 to the solid angle
    
end
SD2 = P*41253.0;
disp(['The fraction of ',num2str(conf2*100),'% to the solid angle is :',num2str(P)]);
disp(['Square degree for ',num2str(conf2*100),'% is :',num2str(SD2)]);
figure;
plot(data(:,1),data(:,2),'m+');
hold on;
% plot(theta,phi,'.r');
% hold on;
plot(xv1,yv1,'k');
hold on
plot(xv2,yv2,'b')
legend('data points',num2str(conf1),num2str(conf2));
xlim([0 2*pi])
ylim([-pi/2 pi/2])
%title('95% confidence contour');
%contour3(X,Y,density,[121.9 121.9])

% End