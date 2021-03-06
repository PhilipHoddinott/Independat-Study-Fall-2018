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


numbSamp=5; % set numbSamp

% sigmaA=linspace(.5,.7,120); % set sigma (km) to go over
% TAdistA=linspace(55,70,120); % set TA dist to go over
sigmaA=linspace(0.001,1,200); % set sigma (km) to go over
TAdistA=linspace(.001,120,300); % set TA dist to go over

coe=coeM(:,1:6);

OrbType='circ';  
eccA=sigmaA;
sigmaOne=1;

[rmsCir,rmsMH,rSigHH] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp,sigmaOne,eccA);
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

figure(42069)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('ecc')
xlabel('\tri M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS');
title(tiS)
colorbar

%% new sec
close all
figure(89)
surf(TAdistA(offSet:end),sigmaA(2:end),rSigHH(2:end,offSet:end));
ylabel('sigma')
xlabel('\tri M (deg)')
zlabel('RMS')
%set(gca,'zscale','log')
tiS=sprintf('RMS');
title(tiS)
colorbar
zlim([-2e4 2e4])
drawnow


keyboard
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
     

function [rmsM,rmsMH,rSigHH] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp,sigmaOne,eccA)
    if strcmp(OrbType,'circ')%OrbType=='circ'
        fgid=101;
    else
        fgid=99;
    end
    figure(fgid)

    h = animatedline('Color','r');
    hTot = animatedline('Color','k');
    hCurr= animatedline('Color','b');
    legend('red: timeLeft est ','black: total time est','blue: current time')
    grid on
    xdrw=1:1:length(sigmaA);
    tdrw=[];
    xc=1;
    


    tic
    fprintf('%.3f\n',toc);
    
    tFlightF=2*pi/(2*pi*rp^(3/2) /sqrt(mu));
    tFlightF=1/tFlightF;
    Tt=2*pi/( (mu^2 *(1-coeM(1,2)^2)^(3/2))/ coeM(1,1)^3);
    
    for sigmaC=1:length(eccA)
        e=eccA(sigmaC);
        inc = 30; % deg
        RAAN = 40;% deg
        w = 70;% deg, arg of perapsis
        rp = 7178.1; % km 

        a=rp*(1-e); % get semi major
        ra=a/(1+e); % get appoapsis
        hM=sqrt(a*(1-e^2)*mu); % get momentum
        TAd=[47,107,138]; % TA
        %% Note in hindsight this is not super nessisary
        for i=1:3 % add given TAs
            TA=TAd(i);
            coeM(i,:)=[hM e RAAN inc w TA a];
        end

        coe=coeM(:,1:6);
        coeLp=coe;
        aReal=coeM(1,7); % set aReal
        
        
        sigma=sigmaOne;
        for TAdistC=1:length(TAdistA)%120
            TAdist=TAdistA(TAdistC);
            TAarr=(0:TAdist:2*TAdist)';

            E=2*atan( sqrt(1-coeM(1,2))/sqrt(1+coeM(1,2)) *tan(.5*TAarr*pi/180));
            MAarr= E-coeM(1,2)*sin(E);
            
            tf21=(Tt/(2*pi))*(MAarr(2)-MAarr(1));%*pi/180;
            tf32=(Tt/(2*pi))*(MAarr(3)-MAarr(2));%*pi/180;
            tf31=(Tt/(2*pi))*(MAarr(3)-MAarr(1));%*pi/180;

            TAarr=TAarr*pi/180;
            [r1T, ~] = sv_from_coe(coeLp(1,:),mu);
            [r2T, v2T] = sv_from_coe(coeLp(2,:),mu);
            [r3T, ~] = sv_from_coe(coeLp(2,:),mu);
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
               
              % v2HTrue=-tf32*( 1/(tf21*tf31) + mu/(12*norm(r1T)^3))*r1T+(tf32-tf21)*( (1/(tf21 *tf32)) + mu/(12*norm(r2T)^3))*r2T+ tf21*( 1/(tf32*tf31) + mu/(12*norm(r3T)^3))*r3T;

               [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);

               coe = coe_from_sv(r2p,v2p,mu);
               a(k,1)=coe(7);
               %coeHH=coe_from_sv(r2p,v2HH,mu);
               coeHH=coe_from_sv(r2,v2HH,mu);
               aHH(k,1)=coeHH(7);
              % v2HH
               %v2T
               %keyboard
               %vDiff(k,1)=v2HH-v2T;
               
            end
           aR=ones(length(a),1)*aReal;
           rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
           rmsHH(sigmaC,TAdistC)=sqrt(mean((aHH(:)-aR).^2));
         % keyboard
            %rSigHH(sigmaC,TAdistC)=v2HH-v2T;
            %mean((aHH(:)-aR)
            rSigHH(sigmaC,TAdistC)=mean( aHH);%v2HH-v2T
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

    
        xc=xc+1;
        drawnow
        
    end
    rmsM=rmsP;
    rmsMH=rmsHH;
end