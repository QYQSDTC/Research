%%% calculate the average CorrelationCoefficients between every two pulsar
a=zeros(1,136);
for i = 1:1:100 %% calculate 100 times
    a(i,:) = CorrelationCoefficient()
end
AC=mean(a);