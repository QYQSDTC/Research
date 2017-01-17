% calculate the confidence elliptic from a 2-D scatter plot, using KDE 
% function 'kde2d'and determine the contour level to use 
% Sep 26, 2016. YW 

clear; 

% load data into workspace 
%results=load('/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc6/resultsSKA_snr0123_AnglesWrapped_v1_May26_2016/summary/H1stats_snr1_loc6_omg3.mat');
results=load('/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc9/resultsSKA_snr0123_loc9/summary/H1stats_snr1_loc9_omg3.mat');

data=[results.estSigParams(:,1), results.estSigParams(:,2)];
Ns=length(results.estSigParams(:,1));  % sample size
N68 = floor(Ns*0.68);
N95 = floor(Ns*0.95);
ra=mean(data(:,1));  % mean value of RA
dec=mean(data(:,2));  % mean value of Dec

% calculate the density by kde2d function
[bandwidth,density,X,Y]=kde2d(data);

% plot the contour to have the coordinates of different contour levels
Ncl=100;  % number of contour level
figure;
[C,h]=contour3(X,Y,density,Ncl);
level=zeros(Ncl,1);  % contour level

Np=zeros(Ncl,1);  % number of date points located inside that contour
cont=1;  % counter
xq=data(:,1);
yq=data(:,2);

for i=1:1:44 %Ncl
    
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

for i=1:1:Ncl   
    if Np(i)<N68
        disp(['Index for 68% level is: ', num2str(i), ', Np= ', num2str(Np(i))...
             ', contour level is: ', num2str(level(i))]);
        break;
    end
end

for i=1:1:Ncl   
    if Np(i)<N95
        disp(['Index for 95% level is: ', num2str(i), ', Np= ', num2str(Np(i))...
             ', contour level is: ', num2str(level(i))]);
        break;
    end
end

%contour3(X,Y,density,[121.9 121.9])

% End