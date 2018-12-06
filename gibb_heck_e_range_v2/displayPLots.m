close all; clear all;
%{
load('slider.mat','bestMult','chrSurf','TAdistA','sigmaA','rmsBest','offSet')


for i=1:14
    figure(i)
    surf(TAdistA(offSet:end),sigmaA(2:end),bestMult(2:end,offSet:end,i));
    title(['ecc = ',num2str(cell2mat(chrSurf(i,2)))])
    ylabel('sigma')
    xlabel('\bigtriangleup M (deg)')
    zlabel('RMS')
    set(gca,'zscale','log')
    colorbar
end
%}

for cTER=1:14
    strLd=sprintf('_wkspc_i_%d.mat',cTER);
    load(strLd);
    figure(cTER)
    title(['ecc = ',num2str(e)])
    subplot(2,2,1)
    hold on
surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end),'FaceAlpha',0.5)
%surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end),'FaceAlpha',0.5)
%surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end))%,'FaceAlpha',0.5)
surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end))%,'FaceAlpha',0.5)
%s1.EdgeColor = 'none';
%s2.EdgeColor = 'none';
ylabel('\sigma')
xlabel('\Delta M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
tiS=sprintf('both plots');
title(tiS)
colorbar
grid on
set(gca,'colorscale','log')
view(17,22)
    %{
    
        surf(TAdistA(offSet:end),sigmaA(2:end),rmsCir(2:end,offSet:end));
        title('gibbs')
    ylabel('sigma')
    xlabel('\bigtriangleup M (deg)')
    zlabel('RMS')
    set(gca,'zscale','log')
    colorbar
    
    %}
    %{
    subplot(2,1,2)
        surf(TAdistA(offSet:end),sigmaA(2:end),rmsMH(2:end,offSet:end));
        title('heck gibbs')
    ylabel('sigma')
    xlabel('\bigtriangleup M (deg)')
    zlabel('RMS')
    set(gca,'zscale','log')
    colorbar
    %}
    subplot(2,2,2)
    surf(TAdistA(offSet:end),sigmaA(2:end),rmsBest(2:end,offSet:end));
    title('best')
    ylabel('\sigma')
    xlabel('\Delta M (deg)')
    zlabel('RMS')
    set(gca,'zscale','log')
    %ylim([0,2])
    %xlim([0 60])
    
    grid on
    colorbar
    set(gca,'colorscale','log')
    %drawnow
    %rmsCir
    %rmsMH
    %311
    %312
    subplot(2,2,3)
    
    vct=[1:1:22];
    %figure(8)
    contour(rmsBest(:,2:end),vct,'ShowText','on')
    ylabel('\sigma')
    xlabel('\Delta M (deg)')
    %ylim([0,2])
    %xlim([0 60])
    grid on
    for ik=1:length(rmsDiff(:,1))
        for j=4:length(rmsDiff(1,:))
            if rmsDiff(ik,j)<0
                pntBet2(ik)=j;
                break;
            end
        end
    end
    pntChr(i)={pntBet2};
    subplot(2,2,4)
    hold on
    for j=1:i
        plot(cell2mat(pntChr(j)),sigmaA,'-s')
        titChr(j)={['e = ',num2str(eArr(j))]};
    end
    legend(titChr,'location','best')
    grid on
    ylabel('sigma')
    xlabel('\delta M')
    %ylim([0,2])
    %xlim([0 60])
    
    drawnow%lose allclea rall;
    
end