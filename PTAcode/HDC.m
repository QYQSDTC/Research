%% calculate the average CorrelationCoefficients between every two pulsar
%% 2017.3.12 QYQ
%% Main

a=zeros(1,136);% there are 17 pulsars in use so there are 16*17/2 coefficients
b=zeros(1,136);
for i = 1:1:2 % calculate 100 times
    [CE,thetaC]=CorrelationCoefficient();
    a(i,:) = CE;
    b(i,:)=thetaC;
end
AC=mean(a);
AT=mean(b);