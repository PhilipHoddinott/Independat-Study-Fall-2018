clear all; close all;
%% Do my best to remove for loops to make vectoried
load('coe_elp','coeM')
mu  = 398600; % mu for earth
rng('default') % For reproducibility
s = rng;
eArr=.0001:.05:.70001
for i=1:length(eArr)
    %% Create orbit
    inc = 30; % deg
    RAAN = 40;% deg
    %e = .1; % ecc for now, e = 0, e = 1.2
    e=eArr(i);
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


    numbSamp=1500; % set numbSamp
    %numbSamp=10; % set numbSamp
    sigmaA=linspace(0,2,100); % set sigma (km) to go over
    TAdistA=linspace(1,60,200); % set TA dist to go over
    coe=coeM(:,1:6);

    OrbType='circ';  
    
    [rmsCir,rmsMH] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp);
    rmsEc=rmsCir;
    offSet=1;
    %% pLots
    close all;

    rmsInd=[];
    rmsBest=[];
    for ik=1:length(rmsCir(:,1))
        for j=1:length(rmsCir(1,:))
            if rmsCir(ik,j)>rmsMH(ik,j)
                rmsBest(ik,j)=rmsMH(ik,j);
                rmsInd(ik,j)=1; % fokr HH
            else
                rmsBest(ik,j)=rmsCir(ik,j);
                rmsInd(ik,j)=0; % for gib
            end
        end
    end
    rmsDiff=rmsCir-rmsMH;

    close all;
    for ik=1:length(rmsDiff(:,1))
        for j=1:length(rmsDiff(1,:))
            if rmsDiff(ik,j)<0
                pntBet2(ik)=j;
                break;
            end
        end
    end
    pntChr(i)={pntBet2};
    figure(12)
    hold on
    for j=1:i
        plot(cell2mat(pntChr(j)),sigmaA,'-s')
        titChr(j)={['e = ',num2str(eArr(j))]};
    end
    legend(titChr)
    grid on
    ylabel('sigma')
    xlabel('delta M')
    pntChr(i)={pntBet2};
    strNm=sprintf('_wkspc_i_%d',i);
    save(strNm)
end

figure(12)
hold on
for j=1:i
    plot(cell2mat(pntChr(j)),sigmaA,'-s')
    titChr(j)={['e = ',num2str(eArr(j))]};
end
legend(titChr)
grid on
ylabel('sigma')
xlabel('delta M')
keyboard
%close all;
%% make gifs
h = figure(2);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for i = 1:length(eArr)
    %strNm=sprintf('wkspc_e_%d',100*eArr(i));
    strNm=sprintf('wkspc_i_%d',i);
    load(strNm,'TAdistA','sigmaA','rmsBest','offSet')

    surf(TAdistA(offSet:end),sigmaA(2:end),rmsBest(2:end,offSet:end));
    ylabel('sigma')
    xlabel('\bigtriangleup M (deg)')
    zlabel('RMS')
    set(gca,'zscale','log')
    %plot([1:5],rmsP)
    tiS=sprintf('RMS best for 3 = %.3f',eArr(i));
    title(tiS)
    colorbar
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
            

function [rmsM,rmsMH] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp)
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