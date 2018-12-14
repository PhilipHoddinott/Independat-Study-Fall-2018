close all; clear all;
offSet=1;
for loopCount=1:3
    str=['wkspcNew_',num2str(loopCount)];
    load(str)
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
    zlabel('RMS-a (km)')
    set(gca,'zscale','log')
    tiS=sprintf('Best RMS-a, ecc = %.3f',e);
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
    zlabel('RMS - a (km)')
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
    
%     
%     figure(5+5*(loopCount-1))
%     iPhilinc=0;
%     hold on
%     
%     surf(MAdistA1(offSet:end),sigmaA1(2:end),rmsCir1(2:end,offSet:end),'FaceColor','interp')%, 'EdgeColor', 'none');
%     ylabel('\sigma (km)')
%     xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
%     zlabel('RMS - a ')
%     set(gca,'zscale','log')
%     tiS=sprintf('Gibbs Method RMS-a, ecc = %.3e',e);
%     %title(tiS)
%     set(gca,'ColorScale','log')
%     colorbar
%     zlim([0.1 10e4])
%     caxis([.1 10e4])
%     lgndcr(1)={'SurfacePlot of Gibbs method'};
%     iPhilinc=iPhilinc+1;
%     
%     myColorMap = jet(length(eArr)); 
%     for iPhil=1:loopCount
%         %line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)), 'Color', 'r', 'LineWidth', 4);
%         line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)),'Color',myColorMap(iPhil,:),'LineWidth', 4);
%         lgndcr(iPhil+iPhilinc)={['e = ',num2str(eArr(iPhil))]};
%     end
%     legend(lgndcr)
%      view(17,22)
%      grid on;
%     drawnow
%     if loopCount==8    
%     %keyboard
%     end
end