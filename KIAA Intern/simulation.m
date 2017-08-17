%% simulate gaussian signal to fit the DM and Residuals
%% QYQ 2017.8.8 in KIAA
%%
clear;
tic;

%% generate pulses
np=-100:0.1:100;
N=length(np);
noise=random('norm',0,1,1,N);
weights=abs(noise);
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
%save('data.txt','xx','yy','-ascii');
% figure
% for i = 1:1:4
%     subplot(2,2,i)
%     plot(t,signal(i,:))
%     hold on
% end
% 
% for i = 2:1:BW
% tao(i)=4150*DM*(ch(i)^(-2)-ch(i-1)^(-2));
% gt=snr*normpdf(np,sum(tao(1:i,:)),0.1);
% g(:,i)=gt;
% signal(:,i)=signal(:,i-1)+gt';
% end
% 
% figure 
% plot(np,signal(:,BW))

% fit_np=np';% to import to origin
% save('fit_x.txt','fit_np','-ascii')
% save('fit_y.txt','signal','-ascii')
% save('weights.txt','weights','-ascii')

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
    options=fitoptions('Weights',weights);
    %ft=fittype('gaussfit(x,a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5,a6,b6,c6,a7,b7,c7,a8,b8,c8,a9,b9,c9,a10,b10,c10,a11,b11,c11,a12,b12,c12,a13,b13,c13,a14,b14,c14,a15,b15,c15,a16,b16,c16,a17,b17,c17,a18,b18,c18,a19,b19,c19,a20,b20,c20,a21,b21,c21,a22,b22,c22,a23,b23,c23,a24,b24,c24,a25,b25,c25,a26,b26,c26,a27,b27,c27,a28,b28,c28,a29,b29,c29,a30,b30,c30,a31,b31,c31)');
    ft=fittype('gauss8');
    y=signal(:,i);
    [p,~]=fit(np(:),y(:),ft,options);
    mu(i)=p.b1;
    plot(p,np,y)
end

ft2=fittype('a*x^(-b)-207.5+84.7189');
f=ch;
t1=mu-tao(BW);
[p,g]=fit(f(:),t1(:),ft2);
error=confint(p);
DM_error=(error(2,1)-error(1,1))/4150;
beta_error=error(2,2)-error(1,2);

p1=predint(p,f,0.95,'functional','on');
figure
plot(p,f,t1)
hold on
plot(f,p1,'m--')
legend('data','fitted','prediction bounds');



%% analytical 
b=p.b;
A=sqrt(fh^(2*b)*fl^(2*b)*(fh^(2*b)*fl-fh*fl^(2*b))*(-1+2*b)^3/((fh^(2*b)*fl-fh*fl^(2*b))^2-fh^(2*b+1)...
    *fl^(2*b+1)*(1-2*b)^2*(log(fh/fl))^2));
dbeta=2^(1/2)*0.1/(4150*DM)*30^(1/2)*snr^(-1)*A;

%save('error.txt','DM_error','beta_error','dbeta','-ascii')
fid=fopen('conclusion.txt','w');
fprintf(fid,'Input Parameters?\n');
fprintf(fid,'DM 500\n beta 2\n BW 100MHz-130MHz\n Channel 31\n');
fprintf(fid,'fit error:\n');
fprintf(fid,'DM error is %8.5f\n',DM_error);
fprintf(fid,'beta error for fit is %8.5f\n',beta_error);
fprintf(fid,'theoretical error for beta is %8.5f',dbeta);
fclose(fid);