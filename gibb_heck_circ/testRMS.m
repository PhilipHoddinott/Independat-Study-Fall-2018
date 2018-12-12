clear all; close all
sigma=1;
numbSamp=100;
a = normrnd(0,sigma,[numbSamp,1]);
a=a+50;
aR=50*ones(numbSamp,1);
rmsA=sqrt(mean((a(:)-aR).^2))