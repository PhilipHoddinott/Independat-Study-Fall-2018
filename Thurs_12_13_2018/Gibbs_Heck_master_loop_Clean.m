%% Gibbs and Herrick-Gibbs Main Code
% By Philip Hoddinott
%% Setup
clear all; close all;
%% Load variables
load('coe_elp','coeM')
mu  = 398600; % mu for earth
rng('default') % For reproducibility
s = rng;


loopCount=1;
eArr=0.001:.001:.005;
load('wkspcNew_1.mat','MAdistA','sigmaA','rmsCir');
   
MAdistA1=MAdistA;
sigmaA1=sigmaA;
rmsCir1=rmsCir;
    %length(eArr)
    %keyboard

for loopCount=1:length(eArr)
    e=eArr(loopCount);
    %% Create orbit
    inc = 30; % deg
    RAAN = 40;% deg
    %e = 0.3; % ecc for now, e = 0, e = 1.2
    %e=0.6;
    %e=0.1;
    %e=0;
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


    numbSamp=5000; % set numbSamp
    TAdistA=linspace(1,30,30*6); % set TA dist to go over
    sigmaA=linspace(0,2,100); % set sigma (km) to go over
    
    OrbType='circ';  
    coe=coeM(:,1:6);

    
    MAdistA=TAdistA;

    for i=1:length(MAdistA)
        Earr(i)=kepler_E(e,MAdistA(i)*pi/180);
    end
    TAdistAReal=2*atan(sqrt((1+e)/(1-e))*tan(Earr/2));
    TAdistA=TAdistAReal*180/pi;

    [rmsCir,rmsMH,rmsGv2,rmsHGv2] = RMS_G_HG_Fast(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp,MAdistA);
    close all;
    str=['wkspcN3_',num2str(loopCount)];
    save(str)
    
    rmsEc=rmsCir;
    offSet=1;
    %% pLots

    figure(1+5*(loopCount-1))
    surf(MAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end),'FaceColor','interp')%, 'EdgeColor', 'none');
    ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS - a ')
    set(gca,'zscale','log')
    tiS=sprintf('Gibbs Method RMS-a, ecc = %.3e',e);
    title(tiS)
    set(gca,'ColorScale','log')
    colorbar
    zlim([0.1 10e4])
    caxis([.1 10e4])
    %strF=[tiS,'.fig'];
    %savefig(1+5*(loopCount-1),strF)
    


    figure(2+5*(loopCount-1))
    surf(MAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end),'FaceColor','interp')%, 'EdgeColor', 'none')
    ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS-a')
    set(gca,'zscale','log')
    tiS=sprintf('Herrick Gibbs RMS-a, ecc = %.3e',e);
    title(tiS)
    set(gca,'ColorScale','log')
    colorbar
    %strF=[tiS,'.fig'];
    %savefig(2+5*(loopCount-1),strF)



    rmsInd=[];
    rmsBest=[];
    for i=1:length(rmsCir(:,1))
        for j=1:length(rmsCir(1,:))
            if rmsCir(i,j)>rmsMH(i,j)
                rmsBest(i,j)=rmsMH(i,j);
                rmsInd(i,j)=1; % for HH
                %rmsBestV2(i,j)=rmsHGv2(i,j);
            else
                rmsBest(i,j)=rmsCir(i,j);
                %rmsBestV2(i,j)=rmsGv2(i,j);
                rmsInd(i,j)=0; % for gib
            end
        end
    end
    rmsDiff=rmsCir-rmsMH;

    z1=rmsCir(2:end,offSet:end);
    z2=rmsMH(2:end,offSet:end);
    x=MAdistA(offSet:end);
    y=sigmaA(2:end);

    zdiff = z1 - z2;
    C = contours(x, y, zdiff, [0 0]);
    
    xL = C(1, 2:end);
    yL = C(2, 2:end);

    zL = interp2(x, y, z1, xL, yL);
    


    figure(3+5*(loopCount-1))
    surf(MAdistA(offSet:end),sigmaA(2:end),rmsBest(2:end,offSet:end),'FaceColor','interp')%, 'EdgeColor', 'none');
    ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS-a')
    set(gca,'zscale','log')
    tiS=sprintf('Best RMS-a, ecc = %.3e',e);
    title(tiS)
    set(gca,'ColorScale','log')
    colorbar
    line(xL, yL, zL, 'Color', 'r', 'LineWidth', 4);
    %strF=[tiS,'.fig'];
    %savefig(3+5*(loopCount-1),strF)

    figure(4+5*(loopCount-1))

    surface(MAdistA(offSet:end),sigmaA(2:end),z1,'FaceColor','interp')%, 'EdgeColor', 'none');

    set(gca,'ColorScale','log')
    colorbar
    zlim([0.1 10e4])
    caxis([.1 10e4])
    hold on

    surf(MAdistA(offSet:end),sigmaA(2:end),z2,'FaceColor','interp')%, 'EdgeColor', 'none');


    ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS - a ')
    set(gca,'zscale','log')
    tiS=sprintf('Both Methods RMS-a, ecc = %.3f',e);
    title(tiS)
    %colorbar
    grid on
    %strF=[tiS,'.fig'];
    %savefig(4+5*(loopCount-1),strF)

    view(17,22)


    line(xL, yL, zL, 'Color', 'r', 'LineWidth', 4);
    
    xLM{loopCount}=xL;
    yLM{loopCount}=yL;
    zLM{loopCount}=zL;
    
    
    figure(5+5*(loopCount-1))
    iPhilinc=0;
    hold on
    
    surf(MAdistA1(offSet:end),sigmaA1(2:end),rmsCir1(2:end,offSet:end),'FaceColor','interp')%, 'EdgeColor', 'none');
    ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS - a ')
    set(gca,'zscale','log')
    tiS=sprintf('Gibbs Method RMS-a, ecc = %.3e',e);
    %title(tiS)
    set(gca,'ColorScale','log')
    colorbar
    zlim([0.1 10e4])
    caxis([.1 10e4])
    lgndcr(1)={'SurfacePlot of Gibbs method'};
    iPhilinc=iPhilinc+1;
    
    myColorMap = jet(length(eArr)); 
    for iPhil=1:loopCount
        %line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)), 'Color', 'r', 'LineWidth', 4);
        line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)),'Color',myColorMap(iPhil,:),'LineWidth', 4);
        lgndcr(iPhil+iPhilinc)={['e = ',num2str(eArr(iPhil))]};
    end
    legend(lgndcr)
     view(17,22)
     grid on;
    drawnow
    if loopCount==8    
    %keyboard
    end
    
end


function [rmsM,rmsMH,rmsGv2,rmsHGv2] = RMS_G_HG_Fast(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp,MAdistA)
    rmsGv2=0;
    rmsHGv2=0;
    rmsP=zeros(length(sigmaA),length(TAdistA));
    rmsHH=zeros(length(sigmaA),length(TAdistA));
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

    coe=coeM(:,1:6);
    coeLp=coe;
    aReal=coeM(1,7); % set aReal
    
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
            E=2*atan( sqrt(1-coeM(1,2))/sqrt(1+coeM(1,2)) *tan(.5*TAarr*pi/180));
            MAarr= E-coeM(1,2)*sin(E);
            
            MAC=MAdistA(TAdistC);
            MACarr=(0:MAC:2*MAC)';
            
            MAarr(3)=MACarr(3)*pi/180;
                       
            tf21=(Tt/(2*pi))*(MAarr(2)-MAarr(1));%*pi/180;
            tf32=(Tt/(2*pi))*(MAarr(3)-MAarr(2));%*pi/180;
            tf31=(Tt/(2*pi))*(MAarr(3)-MAarr(1));%*pi/180;
            %coeM(i,:)=[h e RAAN inc w TA a];
            %taF21=sqrt(coeM(1,7)^3 /mu)*(E(2)-E(1)-coeM(1,2)*sin(E(2)-E(1)));
            for kCount=1:3
                t_hGibbs(kCount)=sqrt(coeM(1,7)^3 /mu)*(E(kCount)-coeM(1,2)*sin(E(kCount)));
            end
            %keyboard

            TAarr=TAarr*pi/180;
            
            

            coeLp(:,6)=TAarr;
            for i=1:3
                [r, v] = sv_from_coe(coeLp(i,:),mu);
                rn = normrnd(0,sigma,[numbSamp,3]);
                rRand(1:numbSamp,1:3)=r(1:3)+rn(1:numbSamp,1:3);%/1000;
                 %rRand(1:numbSamp,1:3)=r(1:3)+0.05*ones(numbSamp,3);%rn(1:numbSamp,1:3);%/1000;
                rMast(:,:,i)=rRand;
            end
            for k=1:numbSamp
               r1=rMast(k,:,1);
               r2=rMast(k,:,2);
               r3=rMast(k,:,3);
               v2HH=-tf32*( 1/(tf21*tf31) + mu/(12*norm(r1)^3))*r1+(tf32-tf21)*( (1/(tf21 *tf32)) + mu/(12*norm(r2)^3))*r2+ tf21*( 1/(tf32*tf31) + mu/(12*norm(r3)^3))*r3;
               
               [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);
                %v2Gm(k,1)=norm(v2p);
               coe = coe_from_sv(r2p,v2p,mu);
               a(k,1)=coe(7);
               coeHH=coe_from_sv(r2,v2HH,mu);
               aHH(k,1)=coeHH(7);
               %v2HM(k,1)=norm(v2HH);
            end
            [~, v2] = sv_from_coe(coeLp(2,:),mu);
           aR=ones(length(a),1)*aReal;
           v2R=ones(length(a),1)*norm(v2);
           
           rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
           rmsHH(sigmaC,TAdistC)=sqrt(mean((aHH(:)-aR).^2));
          if TAdist>38.2*pi/180&&TAdist<38.6*pi/180
              % keyboard;
           end
           if TAdistC==127
              % keyboard;
           end
           
           if rmsP(sigmaC,TAdistC)-rmsHH(sigmaC,TAdistC)<0
               rmsHH(sigmaC,TAdistC+1:length(TAdistA))=10^4;
               rmsP(sigmaC,TAdistC+1:length(TAdistA))=10;
               break;
           end
           
           
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


function [rmsM,rmsMH,rmsGv2,rmsHGv2] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp,MAdistA)
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

    coe=coeM(:,1:6);
    coeLp=coe;
    aReal=coeM(1,7); % set aReal
    
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
            E=2*atan( sqrt(1-coeM(1,2))/sqrt(1+coeM(1,2)) *tan(.5*TAarr*pi/180));
            MAarr= E-coeM(1,2)*sin(E);
            
            MAC=MAdistA(TAdistC);
            MACarr=(0:MAC:2*MAC)';
            
            MAarr(3)=MACarr(3)*pi/180;
                       
            tf21=(Tt/(2*pi))*(MAarr(2)-MAarr(1));%*pi/180;
            tf32=(Tt/(2*pi))*(MAarr(3)-MAarr(2));%*pi/180;
            tf31=(Tt/(2*pi))*(MAarr(3)-MAarr(1));%*pi/180;
            %coeM(i,:)=[h e RAAN inc w TA a];
            %taF21=sqrt(coeM(1,7)^3 /mu)*(E(2)-E(1)-coeM(1,2)*sin(E(2)-E(1)));
            for kCount=1:3
                t_hGibbs(kCount)=sqrt(coeM(1,7)^3 /mu)*(E(kCount)-coeM(1,2)*sin(E(kCount)));
            end
            %keyboard

            TAarr=TAarr*pi/180;
            
            

            coeLp(:,6)=TAarr;
            for i=1:3
                [r, v] = sv_from_coe(coeLp(i,:),mu);
                rn = normrnd(0,sigma,[numbSamp,3]);
                rRand(1:numbSamp,1:3)=r(1:3)+rn(1:numbSamp,1:3);%/1000;
                 %rRand(1:numbSamp,1:3)=r(1:3)+0.05*ones(numbSamp,3);%rn(1:numbSamp,1:3);%/1000;
                rMast(:,:,i)=rRand;
            end
            for k=1:numbSamp
               r1=rMast(k,:,1);
               r2=rMast(k,:,2);
               r3=rMast(k,:,3);
               v2HH=-tf32*( 1/(tf21*tf31) + mu/(12*norm(r1)^3))*r1+(tf32-tf21)*( (1/(tf21 *tf32)) + mu/(12*norm(r2)^3))*r2+ tf21*( 1/(tf32*tf31) + mu/(12*norm(r3)^3))*r3;
               
               [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);
                v2Gm(k,1)=norm(v2p);
               coe = coe_from_sv(r2p,v2p,mu);
               a(k,1)=coe(7);
               coeHH=coe_from_sv(r2,v2HH,mu);
               aHH(k,1)=coeHH(7);
               v2HM(k,1)=norm(v2HH);
            end
            [~, v2] = sv_from_coe(coeLp(2,:),mu);
           aR=ones(length(a),1)*aReal;
           v2R=ones(length(a),1)*norm(v2);
           rmsGv2(sigmaC,TAdistC)=sqrt(mean((v2Gm(:)-v2R).^2));
           rmsHGv2(sigmaC,TAdistC)=sqrt(mean((v2HM(:)-v2R).^2));
           rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
           rmsHH(sigmaC,TAdistC)=sqrt(mean((aHH(:)-aR).^2));
          if TAdist>38.2*pi/180&&TAdist<38.6*pi/180
              % keyboard;
           end
           if TAdistC==127
              % keyboard;
           end
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


function [rmsGv2,rmsHGv2,rmsME] = checkDisc(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp)
    if strcmp(OrbType,'circ')%OrbType=='circ'
        fgid=101;
    else
        fgid=99;
    end
    figure(fgid)
    close(fgid)
    figure(fgid)

    h = animatedline('Color','r');
    hTot = animatedline('Color','k');
    hCurr= animatedline('Color','b');
    legend('red: timeLeft est ','black: total time est','blue: current time')
    grid on
    xdrw=1:1:length(sigmaA);
    tdrw=[];
    xc=1;

    coe=coeM(:,1:6);
    coeLp=coe;
    aReal=coeM(1,7); % set aReal
    
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
            E=2*atan( sqrt(1-coeM(1,2))/sqrt(1+coeM(1,2)) *tan(.5*TAarr*pi/180));
            MAarr= E-coeM(1,2)*sin(E);
                       
            tf21=(Tt/(2*pi))*(MAarr(2)-MAarr(1));%*pi/180;
            tf32=(Tt/(2*pi))*(MAarr(3)-MAarr(2));%*pi/180;
            tf31=(Tt/(2*pi))*(MAarr(3)-MAarr(1));%*pi/180;
            %coeM(i,:)=[h e RAAN inc w TA a];
            %taF21=sqrt(coeM(1,7)^3 /mu)*(E(2)-E(1)-coeM(1,2)*sin(E(2)-E(1)));
            for kCount=1:3
                t_hGibbs(kCount)=sqrt(coeM(1,7)^3 /mu)*(E(kCount)-coeM(1,2)*sin(E(kCount)));
            end
            %keyboard

            TAarr=TAarr*pi/180;
            
            

            coeLp(:,6)=TAarr;
            for i=1:3
                [r, v] = sv_from_coe(coeLp(i,:),mu);
                rn = normrnd(0,sigma,[numbSamp,3]);
                %rRand(1:numbSamp,1:3)=r(1:3)+rn(1:numbSamp,1:3);%/1000;
                 rRand(1:numbSamp,1:3)=r(1:3)+0.05*ones(numbSamp,3);%rn(1:numbSamp,1:3);%/1000;
                rMast(:,:,i)=rRand;
            end
            for k=1:numbSamp
               r1=rMast(k,:,1);
               r2=rMast(k,:,2);
               r3=rMast(k,:,3);
               v2HH=-tf32*( 1/(tf21*tf31) + mu/(12*norm(r1)^3))*r1+(tf32-tf21)*( (1/(tf21 *tf32)) + mu/(12*norm(r2)^3))*r2+ tf21*( 1/(tf32*tf31) + mu/(12*norm(r3)^3))*r3;
               
               [r2p,v2p] = gibbs_Fun(r1,r2,r3,mu);
                v2Gm(k,1)=norm(v2p);
               coe = coe_from_sv(r2p,v2p,mu);
               a(k,1)=coe(7);
               coeHH=coe_from_sv(r2,v2HH,mu);
               aHH(k,1)=coeHH(7);
               v2HM(k,1)=norm(v2HH);
            end
            [~, v2] = sv_from_coe(coeLp(2,:),mu);
           aR=ones(length(a),1)*aReal;
           v2R=ones(length(a),1)*norm(v2);
           
           rmsGv2(sigmaC,TAdistC)=mean((v2Gm(:)-v2R));
           rmsHGv2(sigmaC,TAdistC)=mean((v2HM(:)-v2R));
           rmsME(sigmaC,TAdistC)=MAarr(2)*180/pi;
           rmsP(sigmaC,TAdistC)=sqrt(mean((a(:)-aR).^2));
           rmsHH(sigmaC,TAdistC)=sqrt(mean((aHH(:)-aR).^2));
           %if TAdist==38.356783919597990
           if TAdist>38.2*pi/180&&TAdist<38.6*pi/180
               keyboard;
           end
           if TAdistC==127
               keyboard;
           end
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