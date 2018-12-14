close all; clear all;
load('wkspc_FINAL2')
iPhilinc=0;
%figure(5+5*(loopCount-1))
figure
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
    zlim([0.1 10e3])
    caxis([.1 10e3])
    xlim([0 20])
    lgndcr(1)={'SurfacePlot of Gibbs method'};
    iPhilinc=iPhilinc+1;
 myColorMap = jet(length(eArr)); 
    for iPhil=1:loopCount
        %line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)), 'Color', 'r', 'LineWidth', 4);
        line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)),'Color',myColorMap(iPhil,:),'LineWidth', 8);
        lgndcr(iPhil+iPhilinc)={['e = ',num2str(eArr(iPhil))]};
    end
    legend(lgndcr)
     view(17,22)
     grid on;
     iPhilinc=0;
     figure(2)
     for iPhil=1:loopCount
        %line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)), 'Color', 'r', 'LineWidth', 4);
        line(cell2mat(xLM(iPhil)), cell2mat(yLM(iPhil)), cell2mat(zLM(iPhil)),'Color',myColorMap(iPhil,:),'LineWidth', 2);
        lgndcr(iPhil+iPhilinc)={['e = ',num2str(eArr(iPhil))]};
    end
    legend(lgndcr)
     ylabel('\sigma (km)')
    xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')

    grid on
    close all;
    iPhilinc=0;
    for iPda=1:3
        str=['wkspcN_FINAL_',num2str(iPda)];
        load(str)
        sigArr=[25,50,75,100];%,125,150,175,200
        %surf(MAdistA1(offSet:end),sigmaA1(2:end),rmsCir1(2:end,offSet:end),'FaceColor','interp')%, 'EdgeColor', 'none');
        figure(iPda)
        hold on;
        iPhilinc=0;
        for ik=3:length(sigmaA)%7%length(sigArr)
            %plot(MAdistA(offSet:end),rmsCir(sigArr(ik),offSet:end))
            plot(MAdistA(offSet:end),rmsMH((ik),offSet:end))
            %semilog(MAdistA(offSet:end),rmsCir(sigArr(ik),offSet:end))
            %lgndcr2(ik+iPhilinc)={['\sigma = ',num2str(sigmaA(sigArr(ik))),'(km)']};
            lgndcr2(ik+iPhilinc-2)={['\sigma = ',num2str(sigmaA((ik))),'(km)']};
        end
        strT=sprintf('ecc = %.4f',e);
        %title(strT)
        grid on
        legend(lgndcr2,'Location','southeast')
        set(gca,'yscale','log')
        xlabel('$\bigtriangleup$ M (deg)','Interpreter','latex')
        ylabel('RMS - a (km)')
        ylim([0 10^3])
    end
   
        %plot(MAdistA