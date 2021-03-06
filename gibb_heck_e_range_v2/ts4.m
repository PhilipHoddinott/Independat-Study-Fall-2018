%function [] = ts4()
%close all; clear all;
load('slider.mat','bestMult','chrSurf','TAdistA','sigmaA','rmsBest','offSet')
% Plot different plots according to slider location.
% S.fh = figure('units','pixels',...
%               'position',[300 300 300 300],...
%               'menubar','none',...
%               'name','slider_plot',...
%               'numbertitle','off',...
%               'resize','off');    
S.fh = figure('units','pixels');%,...
              %'position',[300 300 300 300],...
              %'menubar','none',...
              %'name','slider_plot',...
              %'numbertitle','off',...
              %'resize','off');   
S.x = 0:.01:1;  % For plotting.         
S.ax = axes('unit','pix',...
            'position',[20 80 260 210]);
S.LN = plot(S.x,S.x,'r');        
%titlSt = sprintf('e = %.3f', someValue);
%S.t=
%set(handles.sliderLabel, 'String', titlSt);
%S.LN=surf(TAdistA(offSet:end),sigmaA(2:end),bestMult(2:end,offSet:end,1));
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 10 260 30],...
                 'min',1,'max',14,'val',1,...
                 'sliderstep',[1/14 1/14],...
                 'callback',{@sl_call,S});  
function [] = sl_call(varargin)
%close all;
% Callback for the slider.
[h,S] = varargin{[1,3]}  % calling handle and data structure.
%h
a=floor(get(h,'value'))
%(S.sl.Value)
set(S.LN,'ydata',S.x.^get(h,'value'))

%set(S.ax ,'title',sprintf('sld=%.3f',S.sl.Value))
%set(S.ax ,'title',['sld=',num2str(S.sl.Value)])
title(['sld=',num2str(get(h,'value'))])
end
