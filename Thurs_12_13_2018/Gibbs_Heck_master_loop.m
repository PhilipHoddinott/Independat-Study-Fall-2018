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
eArr=0.005*3:.005:.2;
load('wkspcNew_1.mat','MAdistA','sigmaA','rmsCir');
   
MAdistA1=MAdistA;
sigmaA1=sigmaA;
rmsCir1=rmsCir;
    

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


    numbSamp=50; % set numbSamp
    %sigmaA=linspace(0,.5,40); % set sigma (km) to go over
    %TAdistA=linspace(35,42,1600); % set TA dist to go over
    TAdistA=linspace(1,35,35*3); % set TA dist to go over

    %sigmaA=linspace(0,2,400); % set sigma (km) to go over
    sigmaA=linspace(0,2,50); % set sigma (km) to go over
    %TAdistA=linspace(1,60,800); % set TA dist to go over
    OrbType='circ';  
    coe=coeM(:,1:6);

    %EdistA=atan(
    MAdistA=TAdistA;

    for i=1:length(MAdistA)
        Earr(i)=kepler_E(e,MAdistA(i)*pi/180);
    end
    TAdistAReal=2*atan(sqrt((1+e)/(1-e))*tan(Earr/2));
    TAdistA=TAdistAReal*180/pi;
    %
    %{
    OrbType='circ';  
    [rmsGv2,rmsHGv2,rmsME] = checkDisc(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp);

    %% Plot
    close all;
    offSet=1;
    figure(1)
    surf(TAdistA(offSet:end),sigmaA(2:end),rmsGv2(2:end,offSet:end));
    ylabel('\sigma')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS')
    %set(gca,'zscale','log')
    tiS=sprintf('RMS');
    title(tiS)
    colorbar

    figure(2)
    surf(TAdistA(offSet:end),sigmaA(2:end),rmsHGv2(2:end,offSet:end));
    ylabel('\sigma')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS')
    %set(gca,'zscale','log')
    tiS=sprintf('RMS');
    title(tiS)
    colorbar

    figure(3)
    surf(TAdistA(offSet:end),sigmaA(2:end),rmsME(2:end,offSet:end));
    ylabel('\sigma')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS')
    %set(gca,'zscale','log')
    tiS=sprintf('RMS');
    title(tiS)
    colorbar

    keyboard
    %}
    [rmsCir,rmsMH,rmsGv2,rmsHGv2] = RMS_G_HG(numbSamp,sigmaA,TAdistA,coeM,mu,OrbType,rp,MAdistA);
    str=['wkspcNew_short_',num2str(loopCount)];
    save(str)
    %load('wkspcFnl.mat')
    rmsEc=rmsCir;
    offSet=1;
    %% pLots
    %close all;
    %close(101);

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



    rmsInd=[];
    rmsBest=[];
    for i=1:length(rmsCir(:,1))
        for j=1:length(rmsCir(1,:))
            if rmsCir(i,j)>rmsMH(i,j)
                rmsBest(i,j)=rmsMH(i,j);
                rmsInd(i,j)=1; % for HH
                rmsBestV2(i,j)=rmsHGv2(i,j);
            else
                rmsBest(i,j)=rmsCir(i,j);
                rmsBestV2(i,j)=rmsGv2(i,j);
                rmsInd(i,j)=0; % for gib
            end
        end
    end
    rmsDiff=rmsCir-rmsMH;

    z1=rmsCir(2:end,offSet:end);
    z2=rmsMH(2:end,offSet:end);
    x=MAdistA(offSet:end);
    y=sigmaA(2:end);
    %axis vis3d
    % Take the difference between the two surface heights and find the contour
    % where that surface is zero.
    zdiff = z1 - z2;
    C = contours(x, y, zdiff, [0 0]);
    % Extract the x- and y-locations from the contour matrix C.
    xL = C(1, 2:end);
    yL = C(2, 2:end);
    % Interpolate on the first surface to find z-locations for the intersection
    % line.
    zL = interp2(x, y, z1, xL, yL);
    % Visualize the line.


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


    figure(4+5*(loopCount-1))
    %hold on;
    % Define the input grid
    %[x, y] = meshgrid(linspace(-1, 1));
    % Calculate the two surfaces
    %z1 = y.^2 + 2*x;
    %z2 = 2*y.^3 - x.^2;
    % Visualize the two surfaces

    %surface(MAdistA(offSet:end),sigmaA(2:end),z1,'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');

    %surf(MAdistA(offSet:end),sigmaA(2:end),z2, 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
    surface(MAdistA(offSet:end),sigmaA(2:end),z1,'FaceColor','interp')%, 'EdgeColor', 'none');

    set(gca,'ColorScale','log')
    colorbar
    zlim([0.1 10e4])
    caxis([.1 10e4])
    hold on
    %zlim([0.1 10e4])
    %colorbar
    %hold on
    %caxis(([0.1 10e4]));
    surf(MAdistA(offSet:end),sigmaA(2:end),z2,'FaceColor','interp')%, 'EdgeColor', 'none');


    ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
    zlabel('RMS - a ')
    set(gca,'zscale','log')
    tiS=sprintf('Both Methods RMS-a, ecc = %.3f',e);
    title(tiS)
    %colorbar
    grid on

    %surface(x, y, z1, 'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');
    %surface(x, y, z2, 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
    %view(3); 
    %camlight; 
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
    
    myColorMap = jet(8); 
    for iPhil=1:loopCount
        %line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)), 'Color', 'r', 'LineWidth', 4);
        line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)),'Color',myColorMap(iPhil,:),'LineWidth', 4);
        lgndcr(iPhil+iPhilinc)={['e = ',num2str(eArr(iPhil))]};
    end
    legend(lgndcr)
     view(17,22)
    drawnow
    if loopCount==8    
    keyboard
    end
    
end
%save('wksipLp3')

keyboard;
figure(12)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ TA (deg)','Interpreter','latex')
zlabel('RMS-a')
set(gca,'zscale','log')
tiS=sprintf('Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar

mAarr=[];
tf21A=[];
tf32A=[];
tf31A=[];
EarrA=[];
Tt=2*pi/( (mu^2 *(1-coeM(1,2)^2)^(3/2))/ coeM(1,1)^3);
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
    mAarr(TAdistC)=MAarr(2)*180/pi;
    tf21A(TAdistC)=tf21;
    tf32A(TAdistC)=tf32;
    tf31A(TAdistC)=tf31;
    EarrA(TAdistC)=E(2);
    mastM(TAdistC,:)=MAarr.*180/pi;
end
figure(13)
surf(mAarr(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ mAarr (deg)','Interpreter','latex')
zlabel('RMS-a')
set(gca,'zscale','log')
tiS=sprintf(' mAarr Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar

figure(14)
surf(tf21A(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ tf21A (sec)','Interpreter','latex')
zlabel('RMS-a')
set(gca,'zscale','log')
tiS=sprintf('tf21A Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar

figure(15)
surf(tf32A(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ tf32A (deg)','Interpreter','latex')
zlabel('RMS-a')
set(gca,'zscale','log')
tiS=sprintf('tf32A Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar

figure(16)
surf(tf31A(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ tf31A (deg)','Interpreter','latex')
zlabel('RMS-a')
set(gca,'zscale','log')
tiS=sprintf('tf31A Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar

figure(17)
surf(EarrA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ EarrA (deg)','Interpreter','latex')
zlabel('RMS-a')
set(gca,'zscale','log')
tiS=sprintf('EarrA Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar

tdiffA=[tf21A;tf32A;tf31A;];
%tf31A(offSet:end)

figure(18)
hold on
%surf(MAdistA(offSet:end),(1:1:3),tdiffA(:,offSet:end))
plot(MAdistA(offSet:end),tdiffA(1,offSet:end))
plot(MAdistA(offSet:end),tdiffA(2,offSet:end))
plot(MAdistA(offSet:end),tdiffA(3,offSet:end))
ylabel('t21, t32, t31')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
%zlabel('tf')
%set(gca,'zscale','log')
tiS=sprintf('tf all Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
legend('t21','t32','t31')
%colorbar

figure(19)
surf(TAdistA(offSet:end),(1:1:3),tdiffA(:,offSet:end))
ylabel('t21, t32, t31')
xlabel('$\bigtriangleup$ TA (deg)','Interpreter','latex')
zlabel('tf')
%set(gca,'zscale','log')
tiS=sprintf('tf all Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
colorbar


figure(20)
hold on
%surf(MAdistA(offSet:end),(1:1:3),tdiffA(:,offSet:end))
plot(MAdistA(offSet:end),mastM(offSet:end,1))
plot(MAdistA(offSet:end),mastM(offSet:end,2))
plot(MAdistA(offSet:end),mastM(offSet:end,3))
ylabel('m1, m2, m3')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
%zlabel('tf')
%set(gca,'zscale','log')
tiS=sprintf('tf all Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
legend('m1','m2','m3')

figure(21)
hold on
%surf(MAdistA(offSet:end),(1:1:3),tdiffA(:,offSet:end))
plot(TAdistA(offSet:end),tdiffA(1,offSet:end))
plot(TAdistA(offSet:end),tdiffA(2,offSet:end))
plot(TAdistA(offSet:end),tdiffA(3,offSet:end))
ylabel('t21, t32, t31')
xlabel('$\bigtriangleup$ TA (deg)','Interpreter','latex')
%zlabel('tf')
%set(gca,'zscale','log')
tiS=sprintf('tf all Herrick Gibbs RMS-a, ecc = %.3f',e);
title(tiS)
legend('t21','t32','t31')
%colorbar

            
%{
figure(11)
surf(MAdistA(offSet:end),sigmaA(2:end),rmsGv2(2:end,offSet:end));
ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
zlabel('RMS-v2')
set(gca,'zscale','log')
tiS=sprintf('RMS-Gibbs');
title(tiS)
colorbar

figure(12)
surf(MAdistA(offSet:end),sigmaA(2:end),rmsHGv2(2:end,offSet:end))
ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
zlabel('RMS-v2')
set(gca,'zscale','log')
tiS=sprintf('RMS-Herrick');
title(tiS)
colorbar
%}

%{
figure(13)
surf(MAdistA(offSet:end),sigmaA(2:end),rmsBestV2(2:end,offSet:end));
ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
zlabel('RMS-v2')
set(gca,'zscale','log')
tiS=sprintf('RMS-Best');
title(tiS)
colorbar
%}
figure(4)
hold on
s1 = surf(MAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end),'FaceAlpha',0.5)
s2 = surf(MAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end),'FaceAlpha',0.5)
s1.EdgeColor = 'none';
s2.EdgeColor = 'none';

ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
zlabel('RMS')
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('both plots');
title(tiS)
%colorbar
set(gca,'zscale','log')
view(17,22)
grid on

%{
figure(5)
surf(MAdistA(offSet:end),sigmaA(2:end),rmsDiff(2:end,offSet:end));
ylabel('\sigma')
xlabel('\bigtriangleup M (deg)')
zlabel('RMS')
zlim([-300, 1000])
%set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('RMS diff (circ-HH)');
title(tiS)
colorbar
%}

figure(7)
hold on
surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end),'FaceAlpha',0.5)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end),'FaceAlpha',0.5)

ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
zlabel('RMS')
set(gca,'zscale','log')

tiS=sprintf('both plots');
title(tiS)
grid on
view(17,22)

vct=[1:1:22];
figure(8)
contour(rmsBest(:,2:end),vct,'ShowText','on')
ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
grid on

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

%% plot pntB
figure(9)
plot(pntBet2,sigmaA,'-s')
ylabel('\sigma')
xlabel('delta M')
grid on


figure(10)
hold on;
% Define the input grid
%[x, y] = meshgrid(linspace(-1, 1));
% Calculate the two surfaces
%z1 = y.^2 + 2*x;
%z2 = 2*y.^3 - x.^2;
% Visualize the two surfaces
z1=rmsCir(2:end,offSet:end);
z2=rmsMH(2:end,offSet:end);
x=MAdistA(offSet:end);
y=sigmaA(2:end);
surface(MAdistA(offSet:end),sigmaA(2:end),z1,'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');
surf(MAdistA(offSet:end),sigmaA(2:end),z2, 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
ylabel('\sigma')
xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
zlabel('RMS - a ')
set(gca,'zscale','log')
tiS=sprintf('Gibbs Method RMS-a, ecc = %.3f',e);
title(tiS)
%colorbar
grid on

%surface(x, y, z1, 'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');
%surface(x, y, z2, 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
%view(3); 
camlight; 
view(17,22)
%axis vis3d
% Take the difference between the two surface heights and find the contour
% where that surface is zero.
zdiff = z1 - z2;
C = contours(x, y, zdiff, [0 0]);
% Extract the x- and y-locations from the contour matrix C.
xL = C(1, 2:end);
yL = C(2, 2:end);
% Interpolate on the first surface to find z-locations for the intersection
% line.
zL = interp2(x, y, z1, xL, yL);
% Visualize the line.
line(xL, yL, zL, 'Color', 'k', 'LineWidth', 3);
            
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
%{
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
%}

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