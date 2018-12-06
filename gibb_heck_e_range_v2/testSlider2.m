close all; clear all; 
chrSurf={};
bestMult=[];
for i=1:14
    strNm=sprintf('_wkspc_i_%d',i);
    load(strNm,'TAdistA','sigmaA','rmsBest','offSet','e')
    chrSurf(i,1)={rmsBest};
    chrSurf(i,2)={e};
    bestMult(:,:,i)=rmsBest;
end
save('slider.mat','bestMult','chrSurf','TAdistA','sigmaA','rmsBest','offSet')
    
zeta = .5;                           % Damping Ratio
wn = 2;                              % Natural Frequency
sys = tf(wn^2,[1,2*zeta*wn,wn^2]); 

%load('_wkspc_i_14.mat', 'pntChr''TAdistA','sigmaA','offSet')
%f = figure;
%ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
%h = stepplot(ax,sys);
%setoptions(h,'XLim',[0,10],'YLim',[0,2]);
f=figure;
%strNm=sprintf('wkspc_i_%d',i);
%load(strNm,'TAdistA','sigmaA','rmsBest','offSet')
offSet=1;
%h=surf(TAdistA(offSet:end),sigmaA(2:end),rmsBest(2:end,offSet:end));
h=surf(TAdistA(offSet:end),sigmaA(2:end),bestMult(2:end,offSet:end,1));%rmsBest(2:end,offSet:end));
ylabel('sigma')
xlabel('\bigtriangleup M (deg)')
zlabel('RMS')
set(gca,'zscale','log')
%plot([1:5],rmsP)
%tiS=sprintf('RMS best for 3 = %.3f',eArr(i));
%title(tiS)
colorbar


b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',1, 'min',1, 'max',14,'SliderStep',[1, 1]);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','14','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','ecc','BackgroundColor',bgcolor);
            
            
            
b.Callback = @(es,ed) updateSystem(h,surf(TAdistA(offSet:end),sigmaA(2:end),bestMult(2:end,offSet:end,es.Value)));%tf(wn^2,[1,2*(es.Value)*wn,wn^2]))