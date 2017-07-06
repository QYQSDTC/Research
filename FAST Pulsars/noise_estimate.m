% Calculate the jitter noise and radiometer noise for fast pulsars
%% QYQ 2017.6.18
clear;
%% pulsar parameters
fastmsp=load('fastmsp1.mat');
N=length(fastmsp.P);% number of pulsars
p=fastmsp.P*10^(-3);%in seconds
w0=fastmsp.W; % the original data for W in degrees
w=w0.*p/360; %fastmsp.W*10^(-3);%w0.*p/360;% transform w to time scale
s=fastmsp.S*10^(-3);%fastmsp.S*10^(-3);
t= 3600; % the integrate time, here we set it as 1h

%% telescope parameters
Tsf=20;% system Temprature for fast
Tsp=28;% ... for parkes
Gf=16.5; % for fast
Gp=0.8;% for parkes
dff=800*10^6;% the width of frequency for fast
dfp=300*10^6; % ... for parkes

%% jitter noise and radiometer noise
sigmaj=0.2.*w.*sqrt(p/t);% jitter noise
sigmarf=(w*Tsf).*sqrt(w./(p-w))./(Gf*s*sqrt(2*dff*t));% radiometer noise for fast
sigmarp=(w*Tsp).*sqrt(w./(p-w))./(Gp*s*sqrt(2*dfp*t));
sigmaf=sqrt(sigmaj.^2+sigmarf.^2);%total sigma for fast
sigmap=sqrt(sigmaj.^2+sigmarp.^2);% ... for parkes

%% dominant noise
jf=zeros(3,1);%  jitter noise dominant is 1 for fast
jp=zeros(3,1);% ...for parkes

for i=1:1:N
    if sigmaj(i) > sigmarf(i)
        jf(i)=1;
    end
    if sigmaj(i) > sigmarp(i)
        jp(i)=1;
    end
end

%% for fast
p100=0;%number of pulsars with sigma less than 100ns
p200=0;
p500=0;

figure;
subplot(1,2,1);
hold on;
for i=1:1:N
    
    if sigmaf(i) < 5*10^(-7)
        p500=p500+1;
        x=fastmsp.D(i).*cosd(fastmsp.b(i)).*cosd(fastmsp.l(i));
        y=fastmsp.D(i).*cosd(fastmsp.b(i)).*sind(fastmsp.l(i));
        h1=plot(x,y,'.r','markersize',15);
       % legend('500ns');
       % hold on;
    end
    
    
    if sigmaf(i) < 2*10^(-7)
        p200=p200+1;
        x=fastmsp.D(i).*cosd(fastmsp.b(i)).*cosd(fastmsp.l(i));
        y=fastmsp.D(i).*cosd(fastmsp.b(i)).*sind(fastmsp.l(i));
        h2=plot(x,y,'.k','markersize',15);
       % legend('200ns');
       % hold on;
    end
    
    
    if sigmaf(i) < 10^(-7)
        p100=p100+1;
        x=fastmsp.D(i).*cosd(fastmsp.b(i)).*cosd(fastmsp.l(i));
        y=fastmsp.D(i).*cosd(fastmsp.b(i)).*sind(fastmsp.l(i));
        h3=plot(x,y,'.g','markersize',15);
      %  legend('100ns');
      %  hold on;
    end
end
title('FAST Pulsars');
xlabel('GalX');
ylabel('GalY');
legend([h1,h2,h3],'500ns','200ns','100ns');
display(['The number of pulsars in fast less than 100ns is ',num2str(p100)]);
display(['The number of pulsars in fast less than 200ns is ',num2str(p200)]);
display(['The number of pulsars in fast less than 500ns is ',num2str(p500)]);


%% for parkes
p100=0;%number of pulsars with sigma less than 100ns
p200=0;
p500=0;

subplot(1,2,2);
hold on;

for i=1:1:N
    
    if sigmap(i) < 5*10^(-7)
        p500 = p500+1;
        x=fastmsp.D(i).*cosd(fastmsp.b(i)).*cosd(fastmsp.l(i));
        y=fastmsp.D(i).*cosd(fastmsp.b(i)).*sind(fastmsp.l(i));
        h1=plot(x,y,'.r','markersize',15);
       % hold on;
        
    end
    
    
    if sigmap(i) < 2*10^(-7)
        p200=p200+1;
        x=fastmsp.D(i).*cosd(fastmsp.b(i)).*cosd(fastmsp.l(i));
        y=fastmsp.D(i).*cosd(fastmsp.b(i)).*sind(fastmsp.l(i));
        h2=plot(x,y,'.k','markersize',15);
       % hold on;
    end
    
    
    if sigmap(i) < 10^(-7)
        p100=p100+1;
        x=fastmsp.D(i).*cosd(fastmsp.b(i)).*cosd(fastmsp.l(i));
        y=fastmsp.D(i).*cosd(fastmsp.b(i)).*sind(fastmsp.l(i));
        h3=plot(x,y,'.g','markersize',15);
       % hold on;
    end
end


%% figure
title('Parkes Pulsars');
xlabel('GalX');
ylabel('GalY');
legend([h1,h2,h3],'500ns','200ns','100ns');
display(['The number of pulsars in Parkes less than 100ns is ',num2str(p100)]);
display(['The number of pulsars in Parkes less than 200ns is ',num2str(p200)]);
display(['The number of pulsars in Parkes less than 500ns is ',num2str(p500)]);
save('sigma.mat','sigmaj','sigmarf','sigmarp','sigmaf','sigmap');


