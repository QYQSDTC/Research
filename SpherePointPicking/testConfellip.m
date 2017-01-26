%The test file for confellip.m
clear;
results=load('/Users/qianyiqian/matlabprograms/spherepointpicking/H1stats_snr1_loc9_omg3.mat');
XX=results.estSigParams(:,1);
YY=results.estSigParams(:,2);
[SD1,SD2]=confellip(XX,YY,'conf1',0.8,'conf2',0.90);
% %syntax: SD1 is the square degree for conf1
%          SD2 is the square degree for conf2
%          conf1 is confidence level 1
%          conf2 is confidence level 2 
%          if conf1 and conf2 are missing, using the default confidence levels for conf1=0.68, conf2=0.95
