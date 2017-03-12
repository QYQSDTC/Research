%%% calculate the average CorrelationCoefficients between every two pulsar
%% 2017.3.12 QYQ
a=zeros(1,136);
b=zeros(1,136);
for i = 1:1:2 %% calculate 100 times
    [CE,thetaC]=CorrelationCoefficient();
    a(i,:) = CE;
    b(i,:)=thetaC;
end
AC=mean(a);
AT=mean(b);