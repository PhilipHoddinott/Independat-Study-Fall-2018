function [r2,v2] = gibbs_Fun(r1,r2,r3,mu)
%GIBBS_FUN Summary of this function goes here
%   Detailed explanation goes here
deg = pi/180;
%global mu
%% doing the example from the textbook
%load('testr1r2r3.mat')


r1m=norm(r1);
r2m=norm(r2);
r3m=norm(r3);

c12=cross(r1,r2);
c23=cross(r2,r3);
c31=cross(r3,r1);

c23hat=c23/norm(c23);
ur1=(r1/r1m);
%c23hat=[.55667,-.66341,.5]
%ur1=

%chk_r1_c23=dot(,c23hat)


N=r1m*c23 + r2m*c31 + r3m*c12;
Nm=norm(N);
D = c12+c23+c31;
Dm=norm(D);

S=r1*(r2m-r3m)+r2*(r3m-r1m)+r3*(r1m-r2m);

v2=sqrt(mu/(Nm*Dm))*(cross(D,r2)/r2m + S);

end

