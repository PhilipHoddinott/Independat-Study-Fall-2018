clear all; close all;
%% Do my best to remove for loops to make vectoried
load('coe_elp','coeM')
mu  = 398600; % mu for earth
rng('default') % For reproducibility
s = rng;

%% Create orbit
inc = 30; % deg
RAAN = 40;% deg
e = 0;%.7; % ecc for now, e = 0, e = 1.2
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


numbSamp=400; % set numbSamp
%numbSamp=10; % set numbSamp
sigmaA=linspace(0,1,20); % set sigma (km) to go over
TAdistA=linspace(1,60,100); % set TA dist to go over
coe=coeM(:,1:6);

OrbType='circ';  
%[rmsCir] = RMS_COE(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType);
[rmsCir,rmsMH] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp);
rmsEc=rmsCir;
offSet=1;
%% pLots
close all;
figure(1)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end));
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS');
title(tiS)
colorbar

figure(2)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS');
title(tiS)
colorbar


%[rmsCir,rmsMH]
rmsInd=[];
rmsBest=[];
for i=1:length(rmsCir(:,1))
    for j=1:length(rmsCir(1,:))
        if rmsCir(i,j)>rmsMH(i,j)
            rmsBest(i,j)=rmsMH(i,j);
            rmsInd(i,j)=1; % for HH
        else
            rmsBest(i,j)=rmsCir(i,j);
            rmsInd(i,j)=0; % for gib
        end
    end
end
rmsDiff=rmsCir-rmsMH;

figure(3)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsBest(2:end,offSet:end));
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS best');
title(tiS)
colorbar

figure(4)
hold on
s1 = surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end),'FaceAlpha',0.5)
s2 = surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end),'FaceAlpha',0.5)
s1.EdgeColor = 'none';
s2.EdgeColor = 'none';

ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('both plots');
title(tiS)
%colorbar
set(gca,'zscale','log')
view(17,22)


figure(5)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsDiff(2:end,offSet:end));
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
zlim([-300, 1000])
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS diff (circ-HH)');
title(tiS)
colorbar




figure(6)
hold on
[xx yy] = meshgrid(1:10);
colormap([1 0 0;0 0 1]) %red and blue
surf(xx,yy,rand(10),ones(10)); %first color (red)
surf(xx,yy,rand(10)+10,ones(10)+1); %second color(blue)
view(17,22)



figure(7)
hold on
surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end),'FaceAlpha',0.5)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end),'FaceAlpha',0.5)
%s1.EdgeColor = 'none';
%s2.EdgeColor = 'none';
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('both plots');
title(tiS)
view(17,22)

vct=[1:1:22];
figure(8)
contour(rmsBest(:,2:end),vct,'ShowText','on')
ylabel('sigma')
xlabel('\tri M (deg)')
grid on


%pntBet=
for j=1:length(rmsDiff(1,:))
    for i=1:length(rmsDiff(:,1))
        if rmsDiff(i,j)<0
            pntBet(j)=i;
            break;
        end
    end
end


for i=1:length(rmsDiff(:,1))
    for j=1:length(rmsDiff(1,:))
        if rmsDiff(i,j)<0
            pntBet2(i)=j;
            break;
        end
    end
end
     
            
%colorbar
        
%{
OrbType='elp';
[rmsP] = RMS_COE(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType);
rmsEc=rmsP;
%}
%{
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


            
figure(5)
hold on
vp2=[(1:1:9),(10:5:50)];
%contour(TAdistA(offSet:end),sigmaA,rmsP(:,offSet:end),150,'ShowText','on')

contour(TAdistA(offSet:end),sigmaA,rmsCir(:,offSet:end),vp2,'-k','ShowText','on')
%}
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


function [rmsM,rmsMH] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp)
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
    
    tFlightF=2*pi/(2*pi*rp^(3/2) /sqrt(mu));
    tFlightF=1/tFlightF;
    Tt=2*pi/( (mu^2 *(1-coeM(1,2)^2)^(3/2))/ coeM(1,1)^3);
    
    for sigmaC=1:length(sigmaA)
        
        sigma=sigmaA(sigmaC);
        for TAdistC=1:length(TAdistA)%120
            TAdist=TAdistA(TAdistC);
            TAarr=(0:TAdist:2*TAdist)';
            %E=2*atan( sqrt( (1-coeM(1,2))/(1+coeM(1,2))) *tan(TAarr*pi/180));
            E=2*atan( sqrt(1-coeM(1,2))/sqrt(1+coeM(1,2)) *tan(.5*TAarr*pi/180));
            MAarr= E-coeM(1,2)*sin(E);
            
            %MAarr=2*atan( sqrt( (1-coeM(1,2))/(1+coeM(1,2)) *tan(TAarr*pi/180))) -coeM(1,2)*sqrt(1-coeM(1,2)^2)*sin(TAarr*pi/180)./(1+coeM(1,2)*cos(TAarr*pi/180));
            
            tf21=(Tt/(2*pi))*(MAarr(2)-MAarr(1));%*pi/180;
            tf32=(Tt/(2*pi))*(MAarr(3)-MAarr(2));%*pi/180;
            tf31=(Tt/(2*pi))*(MAarr(3)-MAarr(1));%*pi/180;
            %MAarr2=2*atan( sqrt( (1-coeM(1,2))/(1+coeM(1,2)) *tan(TAarr*pi/180))) -coeM(1,2)*sqrt(1-coeM(1,2)^2)*sin(TAarr*pi/180)./(1+coeM(1,2)*cos(TAarr*pi/180));
            
            %tf212=(Tt/(2*pi))*(MAarr2(2)-MAarr2(1));%*pi/180;
            %tf322=(Tt/(2*pi))*(MAarr2(3)-MAarr2(2));%*pi/180;
            %tf312=(Tt/(2*pi))*(MAarr2(3)-MAarr2(1));%*pi/180;
           
            %keyboard
            %pause
            %pause
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
               v2HH=-tf32*( 1/(tf21*tf31) + mu/(12*norm(r1)^3))*r1+(tf32-tf21)*( (1/(tf21 *tf32)) + mu/(12*norm(r2)^3))*r2+ tf21*( 1/(tf32*tf31) + mu/(12*norm(r3)^3))*r3;
               
               %tf21
               %tf32
               %tf31
               %v2HH= -
               [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);
               %v2p
               %pause
               coe = coe_from_sv(r2p,v2p,mu);
               a(k,1)=coe(7);
               coeHH=coe_from_sv(r2p,v2HH,mu);
               aHH(k,1)=coeHH(7);
            end
           aR=ones(length(a),1)*aReal;
           rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
           rmsHH(sigmaC,TAdistC)=sqrt(mean((aHH(:)-aR).^2));
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
    rmsMH=rmsHH;
end