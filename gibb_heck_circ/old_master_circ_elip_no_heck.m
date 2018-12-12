clear all; close all;
%% Do my best to remove for loops to make vectoried
load('coe_elp','coeM')
mu  = 398600; % mu for earth
rng('default') % For reproducibility
s = rng;

%% Create orbit
inc = 30; % deg
RAAN = 40;% deg
e = .4; % ecc for now, e = 0, e = 1.2
w = 70;% deg, arg of perapsis
rp = 7178.1; % km 

a=rp*(1-e); % get semi major
ra=a/(1+e); % get appoapsis
h=sqrt(a*(1-e^2)*mu); % get momentum
TAd=[47,107,138]; % TA
%% Note in hindsight this is not super nessisary
for i=1:3 % add given TAs
    TA=TAd(i);
    coeM(i,:)=[h e RAAN inc w TA a];
end


numbSamp=200; % set numbSamp
%numbSamp=10; % set numbSamp
sigmaA=linspace(0,10,50); % set sigma (km) to go over
TAdistA=linspace(1,120,120); % set TA dist to go over
coe=coeM(:,1:6);
OrbType='elp';

MAdistA=TAdistA;

for i=1:length(MAdistA)
    Earr(i)=kepler_E(e,MAdistA(i)*pi/180);
end
TAdistAReal=2*atan(sqrt((1+e)/(1-e))*tan(Earr/2));
TAdistA=TAdistAReal*180/pi;

[rmsP] = RMS_COE(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType);
rmsEc=rmsP;
%% Plot Section
%save('wksp10')
offSet=20;%20
figure(1)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsP(2:end,offSet:end));
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS');
title(tiS)
colorbar

%x=[20,17,23]; xmax=20; x(x>xmax)=xmax
rmsPMax=100;
rmsP_clip=rmsP;
rmsP_clip(rmsP_clip>rmsPMax) = rmsPMax;
figure(11)
surf(TAdistA,sigmaA(2:end),rmsP_clip(2:end,:))
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS clip for  %d samples',numbSamp);
title(tiS)
colorbar

figure(12)
Z=rmsP(2:end,offSet:end);
x=TAdistA(offSet:end);
y=sigmaA(2:end);
[dfdx,dfdy] = gradient(Z);
surf(x,y,Z,sqrt(dfdx.^2 + dfdy.^2))
colorbar
tiS=sprintf('grad RMS for  %d samples',numbSamp);
title(tiS)

figure(13)
Z=rmsP_clip;
x=TAdistA;
y=sigmaA;
[dfdx,dfdy] = gradient(Z);
surf(x,y,Z,sqrt(dfdx.^2 + dfdy.^2))
colorbar
tiS=sprintf('grad RMS clip for  %d samples',numbSamp);
title(tiS)


figure(14)
subplot(2,1,1);
surf(TAdistA,sigmaA(2:end),rmsP_clip(2:end,:))
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS clip for  %d samples',numbSamp);
title(tiS)
colorbar

subplot(2,1,2);
Z=rmsP_clip;
x=TAdistA;
y=sigmaA;
[dfdx,dfdy] = gradient(Z);
surf(x,y,Z,sqrt(dfdx.^2 + dfdy.^2))
colorbar
tiS=sprintf('grad RMS clip for  %d samples',numbSamp);
title(tiS)

figure(2)
vp1=[(1:1:9),(10:2:28),(30:5:100)];
%contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),150,'ShowText','on')
contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),vp1,'ShowText','on')
%contour3(TAdistA,sigmaA,rmsP)
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMScount');
title(tiS)
vp=[1,10,50,100,500];
offSet=1;
figure(3)
hold on;
contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),vp,'ShowText','on')
%contour3(TAdistA,sigmaA,rmsP)
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
grid on
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('cotours for %d samples',numbSamp);
title(tiS)


figure(4)
hold on;
contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),vp,'-r','ShowText','on')

%contour3(TAdistA,sigmaA,rmsP)
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
grid on
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('cotours for %d samples',numbSamp);
title(tiS)




inc = 30; % deg
RAAN = 40;% deg
e = 0; % ecc for now, e = 0, e = 1.2
w = 70;% deg, arg of perapsis
rp = 7178.1; % km 

a=rp*(1-e);
ra=a/(1+e);

h=sqrt(a*(1-e^2)*mu);
TAd=[20,60,100];
for i=1:3
    TA=TAd(i);
    coeM(i,:)=[h e RAAN inc w TA a];
end
    %coe = [h e RA incl w TA a];
    %coe
    OrbType='circ';
    drawnow;
[rmsCir] = RMS_COE(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType);
%{
coe=coeM(:,1:6);
coeLp=coe;
aReal=coeM(1,7);
for sigmaC=1:length(sigmaA)
    sigma=sigmaA(sigmaC);
    for TAdistC=1:length(TAdistA)%120
        TAdist=TAdistA(TAdistC);
        TAarr=(0:TAdist:2*TAdist)';
        TAarr=TAarr*pi/180;
        %coeLp=coe;
        coeLp(:,6)=TAarr;
        for i=1:3
            [r, v] = sv_from_coe(coeLp(i,:),mu);
            rn = normrnd(0,sigma,[numbSamp,3]);
            rRand(1:numbSamp,1:3)=r(1:3)+rn(1:numbSamp,1:3);%/1000;
            %{
            for k = 1:numbSamp

                for ik=1:3
                    rRand(k,ik)=r(ik)+rn(k,ik);
                end
            end
            %}
            %whos
            %keyboard
            rMast(:,:,i)=rRand;

        end
       for k=1:numbSamp
           r1=rMast(k,:,1);
           r2=rMast(k,:,2);
           r3=rMast(k,:,3);
           [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);
           coe = coe_from_sv(r2p,v2p,mu);
           a(k,1)=coe(7);
       end
       aR=ones(length(a),1)*aReal;
       rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
       %whos
       %keyboard
    

    end
    fprintf('%d of %d v2\n',sigmaC,length(sigmaA));
end
rmsCir=rmsP;
%}
contour(TAdistA(offSet:end),sigmaA,rmsCir(:,offSet:end),vp,'-k','ShowText','on')
%save('wkspaceCirc')

            %{
figure(3)
[X,Y,Z] = peaks;
v = [1,1];
contour(X,Y,Z,v)

figure(4)
v = [1,2];
contour(X,Y,Z,v)
    

%}
figure(5)
hold on
vp2=[(1:1:9),(10:5:50)];
%contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),150,'ShowText','on')
contour(TAdistA(offSet:end),sigmaA,rmsEc(:,offSet:end),vp,'-r','ShowText','on')
contour(TAdistA(offSet:end),sigmaA,rmsCir(:,offSet:end),vp,'-k','ShowText','on')
title('red ecc, black circ')
xlabel('mean anom')
ylabel('\sigma, (km)')
grid on

figure(6)
hold on

%contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),150,'ShowText','on')
contour(TAdistA(offSet:end),sigmaA,rmsEc(:,offSet:end),vp2,'-r','ShowText','on')
contour(TAdistA(offSet:end),sigmaA,rmsCir(:,offSet:end),vp2,'-k','ShowText','on')
title('red ecc, black circ')
xlabel('mean anom')
ylabel('\sigma, (km)')
grid on

function [rmsM] = RMS_COE(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType)
    if strcmp(OrbType,'circ')%OrbType=='circ'
        fgid=101;
    else
        fgid=99;
    end
    figure(fgid)
    
        %addpoints(h,xdrw(xc),tdrw(xc),'-r');
        %addpoints(hTot,xdrw(xc),tTot,'-k');
        %addpoints(hCurr,xdrw(xc),tr,'-b');
    h = animatedline('Color','r');
    hTot = animatedline('Color','k');
    hCurr= animatedline('Color','b');
    legend('red: timeLeft est ','black: total time est','blue: current time')
    grid on
    xdrw=1:1:length(sigmaA);
    tdrw=[];
    xc=1;
    %numbSamp=100; % set numbSamp
    %sigmaA=linspace(0,20,200/2); % set sigma (km) to go over
    %TAdistA=linspace(1,120,120/2); % set TA dist to go over
    coe=coeM(:,1:6);
    coeLp=coe;
    aReal=coeM(1,7); % set aReal
    %fprintf('%.3f\n',tic);
    tic
    fprintf('%.3f\n',toc);
    
    for sigmaC=1:length(sigmaA)
        
        sigma=sigmaA(sigmaC);
        for TAdistC=1:length(TAdistA)%120
            TAdist=TAdistA(TAdistC);
            TAarr=(0:TAdist:2*TAdist)';
            TAarr=TAarr*pi/180;

            coeLp(:,6)=TAarr;
            for i=1:3
                [r, v] = sv_from_coe(coeLp(i,:),mu);
                rn = normrnd(0,sigma,[numbSamp,3]);
                rRand(1:numbSamp,1:3)=r(1:3)+rn(1:numbSamp,1:3);%/1000;
                rMast(:,:,i)=rRand;
            end
            for k=1:numbSamp
               r1=rMast(k,:,1);
               r2=rMast(k,:,2);
               r3=rMast(k,:,3);
               [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);
               coe = coe_from_sv(r2p,v2p,mu);
               a(k,1)=coe(7);
            end
           aR=ones(length(a),1)*aReal;
           rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
        end
        fprintf('%d of %d, for %s,  ',sigmaC,length(sigmaA),OrbType);
        tr=toc;
        pctR=sigmaC/length(sigmaA);
        tTot=tr/pctR;
        tLeft=tTot-tr;
        fprintf('%.1f sec elaps, %.1f est tot, %.1f left\n',tr,tTot,tLeft);
        tdrw=[tdrw,tLeft];
        addpoints(h,xdrw(xc),tdrw(xc));
        addpoints(hTot,xdrw(xc),tTot);
        addpoints(hCurr,xdrw(xc),tr);
          % hTot = animatedline;
    %hCurr= animatedline;
    
        xc=xc+1;
        drawnow
        
    end
    rmsM=rmsP;
end