clear all; close all;
load('coe_elp','coeM')
mu  = 398600; % mu for earth
rng('default') % For reproducibility
s = rng;
%% Create orbit
inc = 30; % deg
RAAN = 40;% deg
%e = .1; % ecc for now, e = 0, e = 1.2
e=0;
%e=eArr(i);
w = 70;% deg, arg of perapsis
rp = 7178.1; % km 

a=rp*(1-e); % get semi major
ra=a/(1+e); % get appoapsis
h=sqrt(a*(1-e^2)*mu); % get momentum
TAd=[47,107,138]; % TA
%% Note in hindsight this is not super nessisary
for ij=1:3 % add given TAs
TA=TAd(ij);
coeM(ij,:)=[h e RAAN inc w TA a];
end



sigmaA=linspace(0,2,100); % set sigma (km) to go over
%TAdistA=linspace(1,50,150); % set TA dist to go over
coe=coeM(:,1:6);

OrbType='circ';  

jd1 = juliandate(2004,10,10,12,21,0);
jd2 = juliandate(2004,10,10,12,21,1);
jd2-jd1;
TAdistA=1;
%[rmsCir,rmsMH] = RMS_G_HG(1,0,TAdistA,coeM,mu,OrbType,rp)

Tt=2*pi/( (mu^2 *(1-coeM(1,2)^2)^(3/2))/ coeM(1,1)^3);
TAdist=14;
TAarr=(0:TAdist:2*TAdist)';
%E=2*atan( sqrt( (1-coeM(1,2))/(1+coeM(1,2))) *tan(TAarr*pi/180));
E=2*atan( sqrt(1-coeM(1,2))/sqrt(1+coeM(1,2)) *tan(.5*TAarr*pi/180));
MAarr= E-coeM(1,2)*sin(E);

%MAarr=2*atan( sqrt( (1-coeM(1,2))/(1+coeM(1,2)) *tan(TAarr*pi/180))) -coeM(1,2)*sqrt(1-coeM(1,2)^2)*sin(TAarr*pi/180)./(1+coeM(1,2)*cos(TAarr*pi/180));

tf21=(Tt/(2*pi))*(MAarr(2)-MAarr(1));%*pi/180;
tf32=(Tt/(2*pi))*(MAarr(3)-MAarr(2));%*pi/180;
tf31=(Tt/(2*pi))*(MAarr(3)-MAarr(1));%*pi/180;

coeLp=coeM(:,1:6);
TAarr=TAarr*pi/180;



coeLp(:,6)=TAarr;
[r1, v] = sv_from_coe(coeLp(1,:),mu);
[r2, v] = sv_from_coe(coeLp(2,:),mu);
[r3, v] = sv_from_coe(coeLp(3,:),mu);
[hgbiisv2,theta,theta1,copa,error] = hgibbs2(r1*1000,r2*1000,r3*1000,tf21,tf31,tf32);
hgbiisv2
v2HH=-tf32*( 1/(tf21*tf31) + mu/(12*norm(r1)^3))*r1+(tf32-tf21)*( (1/(tf21 *tf32)) + mu/(12*norm(r2)^3))*r2+ tf21*( 1/(tf32*tf31) + mu/(12*norm(r3)^3))*r3
[gibV2, ierr] = gibbs_Fun(r1, r2, r3,mu);
%gibV2
