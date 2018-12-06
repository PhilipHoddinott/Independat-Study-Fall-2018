function [] = ts6()
%{
%function scrollplotdemo
%%%%%Generate and plot data
x=0:1e-2:2*pi;
y=sin(x);
dx=2;
% dx is the width of the axis 'window'
a=gca;
p=plot(x,y);
%%%%%Set appropriate axis limits and settings
set(gcf,'doublebuffer','on');
% This avoids flickering when updating the axis
set(a,'xlim',[0 dx]);
set(a,'ylim',[min(y) max(y)]);
%%%%%Generate constants for use in uicontrol initialization
pos=get(a,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
% This will create a slider which is just underneath the axis
% but still leaves room for the axis labels above the slider
xmax=max(x);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
% Setting up callback string to modify XLim of axis (gca)
% based on the position of the slider (gcbo)
%%%%%Creating Uicontrol
h=uicontrol('style','slider',...
      'units','normalized','position',Newpos,...
      'callback',S,'min',0,'max',xmax-dx);
%}
h = surf(peaks(20));
CD = get(h,'cdata');
Min = min(CD(:));
Max = max(CD(:));
cb=colorbar;
DeltaC = (Max-Min)/10;
set(cb,'units','normalized','ylim',Min+[0 1]);
Pos=get(cb,'position');
S=['set(findobj(gcf,''tag'',''Colorbar''),''YLim'',get(gcbo,''Value'')+[0 ' num2str(DeltaC) '])'];
uic=uicontrol('style','slider','units','normalized',...
      'position',[Pos(1)+0.11 Pos(2) 0.04 Pos(4)],...
      'min',Min,'max',Max-DeltaC,'value',Min,...
      'callback',S);
end

