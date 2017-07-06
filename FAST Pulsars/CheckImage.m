% calculate the network response matrix and its condition number
% Ref. CQG 23 (2006) S673-685
% Oct. 2, 2014

clear;

%alphaP=zeros(Np,1);
%deltaP=zeros(Np,1);
%distP=zeros(Np,1);  % (parallax) distance from SSB to pulsars, from mas/pc to ly
%kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,

% retrive the configuration of a PTA
%[PTAconf, pname] = PTA_8();
%[PTAconf, pname] = PTA_17();
for t = 100:100:1000
    filename=strcat('ska-',num2str(t),'.mat');
PTAconf=load(filename);
Np=length(PTAconf.alphaP);  % number of detectors/pulsars
alphaP=PTAconf.alphaP;
deltaP=PTAconf.deltaP;
kp=PTAconf.kp;

A=zeros(Np,2);  % network response matrix: Np pulsars, 2 GW polarizations +,x
Fp=zeros(Np,1);  % antenna pattern func for +,  range vector of A
Fc=zeros(Np,1);  % antenna pattern func for x

Na=200;
Nd=200;
alpha=zeros(Na,1);  % ra of source
delta=zeros(Nd,1);  % dec of source
ks=zeros(1,3);
Acond=zeros(Na,Nd);  % condition number of A at a sky location
Inten=zeros(Na,Nd);  % sqrt(Fp^2+Fc^2)

for j=1:1:Na
    alpha(j)=(j-1)*2*pi/Na;
    
    for k=1:1:Nd
        delta(k)=(k-1)*pi/Nd-pi/2;
        
        ks(1)=cos(delta(k))*cos(alpha(j));
        ks(2)=cos(delta(k))*sin(alpha(j));
        ks(3)=sin(delta(k));
        
        for i=1:1:Np
            alphatilde=alpha(j)-alphaP(i);
            theta=acos(ks*kp(i,:)');
            
            Pp=-cos(deltaP(i))^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta(k))^2)...
                +sin(deltaP(i))^2*cos(delta(k))^2-0.5*sin(2*deltaP(i))*cos(alphatilde)*sin(2*delta(k));
            
            Pc=2*cos(deltaP(i))*sin(alphatilde)*(cos(deltaP(i))*cos(alphatilde)*sin(delta(k))...
                -sin(deltaP(i))*cos(delta(k)));
            
            % + polarization
            Fp(i)=Pp/(1-cos(theta));
            A(i,1)=Fp(i);
            
            % x polarization
            Fc(i)=Pc/(1-cos(theta));
            A(i,2)=Fc(i);
            
            %Inten(j,k)=Inten(j,k) + sqrt(Fp(i)^2 + Fc(i)^2);
                     
        end
        
        Inten(j,k)=sqrt(Fp(15)^2 + Fc(15)^2);
        % calculate the condition number
        Acond(j,k)=cond(A,2);  % cond(A,inf)
        
        
    end
end

[v,ind]=max(Acond);
[v1,ind1]=max(max(Acond));
disp(['The largest element in Acond is: ' num2str(v1) ' at (' num2str(ind(ind1)) ',' num2str(ind1) ')']);
disp(['The values for ra and dec are: ', num2str(alpha(ind(ind1))),' and ',num2str(delta(ind1))]);

% skymap of the condition number
figure
surf(alpha,delta,Acond');
%map=parula(1024);
colormap;
caxis([1.00 1.90]);
title(['Skymap for the condition number of A for ', num2str(Np), ' pulsars in ska']);
%surf(alpha,delta,log10(Acond)');  % log10 plot
shading flat;
view(2);
%gird off;
% figure
% surf(alpha,delta,Inten');
% title(['Skymap for Inten for ', num2str(Np), ' pulsars']);
% % imagesc
end

% End