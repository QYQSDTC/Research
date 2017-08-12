%% simulate gaussian signal to fit the DM and Residuals
%% QYQ 2017.8.8 in KIAA
%%
clear;
tic;

%% generate pulses
np=-100:0.1:100;
N=length(np);
noise=random('norm',0,1,1,N);
fl=100;% MHz frequency down limit 140MHz
ch=100:1:130;% different channel
BW=length(ch);

fh=130;% MHz up limit 180MHz
DM=500;
snr=5;
tao=zeros(BW,1);
g=zeros(N,BW);
signal=zeros(N,BW);
t=zeros(N,BW);
t(:,1)=0:0.1:200;
g0=snr*normpdf(np,0,0.1);
g(:,1)=g0;
signal0=noise+g0;
signal(:,1)=signal0;
%%
for i = 2:1:BW
    t(:,i)=(i-1)*200:0.1:i*200;
tao(i)=4150*DM*(ch(i)^(-2)-ch(1)^(-2));
gt=snr*normpdf(np,tao(i),0.1);
g(:,i)=gt;
%signal(i,:)=signal(1,:)+gt;
signal(:,i)=noise+gt;
end

xx=reshape(t,[1 BW*N]);
yy=reshape(signal,[1 BW*N]);
% figure
% for i = 1:1:4
%     subplot(2,2,i)
%     plot(t,signal(i,:))
%     hold on
% end
% 
% for i = 2:1:4
% tao(i)=4150*DM*(ch(i)^(-2.1)-ch(i-1)^(-2.1));
% gt=snr*normpdf(t-sum(tao(1:i,:)),0,0.1);
% g(i,:)=gt;
% signal(i,:)=signal(i-1,:)+gt;
% end
% 
% figure 
% plot(t,signal(4,:))

%% curve fitting
% ft=fittype('gauss1');
% y=signal(4,:);
% [p,gof]=fit(t(:),y(:),ft);
% % figure
% % plot(p,t,y,'-b')
% 
% ft1=fittype('gauss1');
% y=signal(1,:);
% [p1,gof1]=fit(t(:),y(:),ft1);
% % figure
% % plot(p1,t,y,'-b')
% 
% ft2=fittype('gauss1');
% y=signal(2,:);
% [p2,gof2]=fit(t(:),y(:),ft2);
% 
% ft3=fittype('gauss1');
% y=signal(3,:);
% [p3,gof3]=fit(t(:),y(:),ft3);
mu=zeros(BW,1);
for i=1:1:BW
    %options=fitoptions('Weight',noise);
    ft=fittype('gauss1');
    y=signal(:,i);
    [p,~]=fit(np(:),y(:),ft);
    mu(i)=p.b1;
end

ft2=fittype('a*x^(-b)-207.5+84.7189');
f=ch;
t1=mu-tao(BW);
[p,g]=fit(f(:),t1(:),ft2);
figure
plot(p,f,t1,'-b')


%% analytical 
b=p.b;
A=sqrt(fh^(2*b)*fl^(2*b)*(fh^(2*b)*fl-fh*fl^(2*b))*(-1+2*b)^3/((fh^(2*b)*fl-fh*fl^(2*b))^2-fh^(2*b+1)...
    *fl^(2*b+1)*(1-2*b)^2*(log(fh/fl))^2));
dbeta=2^(1/2)*0.1/(4150*DM)*30^(1/2)*snr^(-1)*A;

