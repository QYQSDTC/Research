% Calculate the confidence elliptic for a 2-D scatter plot, use KDE
% function 'kde2d' and determine the contour level to use
% Sep 26, 2016. YW
% Jan 18, 2017. Yi-Qian Qian added contour and Monte Carlo simulation 

clear;

% load data into workspace
%results=load('/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc6/resultsSKA_snr0123_AnglesWrapped_v1_May26_2016/summary/H1stats_snr1_loc6_omg3.mat');
%results=load('F:\GitProject\MatlabPrograms\SpherePointPicking/H1stats_snr1_loc9_omg3.mat');
results=load('/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc9/resultsSKA_snr0123_loc9/summary/H1stats_snr1_loc9_omg3.mat');

data=[results.estSigParams(:,1), results.estSigParams(:,2)]; % x, y coords 
Ns=length(results.estSigParams(:,1));  % sample size
N68 = floor(Ns*0.68);  % contour level
N95 = floor(Ns*0.95);  % contour level
ra=mean(data(:,1));  % mean value of RA
dec=mean(data(:,2));  % mean value of Dec

% calculate the density by 'kde2d' function
[bandwidth,density,X,Y]=kde2d(data);

% plot the contour to have the coordinates of different contour levels
Ncl=100;  % number of contour level
figure;
[C,h]=contour3(X,Y,density,Ncl);
level=zeros(Ncl,1);  % contour level

Np=zeros(Ncl,1);  % number of date points located inside that contour
Np68=zeros(3,1);  % number of random sphere points located inside that contour
Np95=zeros(3,1);  % number of random sphere points located inside that contour
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

% generate the N68 contour
for i=1:1:Ncl
    if Np(i)<N68
        disp(['Index for 68% level is: ', num2str(i), ', Np= ', num2str(Np(i))...
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
    xv=C(1,int:fin)';
    yv=C(2,int:fin)';
    
    in=inpolygon(theta,phi,xv,yv);
    Np68(i)=sum(in);
    P = Np68(i)/n ; %calculate the fraction of N68 to the solid angle
end
disp(['The fraction of N68 to the solid angle is :',num2str(P)]);
disp(['Square degree for N68 is :',num2str(P*41253.0)]);
figure;
plot(data(:,1),data(:,2),'g+');
hold on;
plot(theta,phi,'.r');
hold on;
plot(xv,yv,'k');
xlim([0 2*pi])
ylim([-pi/2 pi/2])
title('68% confidence contour');

% generate the N95 contour
for i=1:1:Ncl
    if Np(i)<N95
        disp(['Index for 95% level is: ', num2str(i), ', Np= ', num2str(Np(i))...
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
    xv=C(1,int:fin)';
    yv=C(2,int:fin)';
    
    in=inpolygon(theta,phi,xv,yv);
    Np95(i)=sum(in);
    P = Np95(i)/n ; % calculate the fraction of N95 to the solid angle
    
end
disp(['The fraction of N95 to the solid angle is :',num2str(P)]);
disp(['Square degree for N95 is :',num2str(P*41253.0)]);
figure;
plot(data(:,1),data(:,2),'g+');
hold on;
plot(theta,phi,'.r');
hold on;
plot(xv,yv,'k');
xlim([0 2*pi])
ylim([-pi/2 pi/2])
title('95% confidence contour');
%contour3(X,Y,density,[121.9 121.9])

% End