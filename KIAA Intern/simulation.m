%% simulate gaussian signal to fit the DM and Residuals
%% QYQ 2017.8.8 in KIAA
%%
clear;
tic;

%% generate pulses
t=-100:0.1:100;
N=length(t);
noise=random('norm',0,1,1,N);
fl=100;% MHz frequency down limit 140MHz
ch=[100 110 120 130];% different channel
fh=130;% MHz up limit 180MHz
DM=500;
snr=10;
tao=zeros(4,1);
g=zeros(4,N);
signal=zeros(4,N);
g0=snr*normpdf(t,0,0.1);
g(1,:)=g0;
signal0=noise+g0;
signal(1,:)=signal0;

for i = 2:1:4
tao(i)=4150*DM*(ch(i)^(-2.1)-ch(1)^(-2.1));
gt=snr*normpdf(t-tao(i),0,0.1);
g(i,:)=gt;
%signal(i,:)=signal(1,:)+gt;
signal(i,:)=noise+gt;
end

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
ft=fittype('gauss1');
y=signal(4,:);
[p,gof]=fit(t(:),y(:),ft);
figure
plot(p,t,y,'-b')

