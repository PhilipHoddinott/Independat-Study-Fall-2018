%{

function testSlider3
x = 1:10;
hplot = plot(x,0*x);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,x,hplot));
end
function makeplot(hObject,event,x,hplot)
n = get(hObject,'Value');
set(hplot,'ydata',x.^n);
drawnow;
end
%}
function testSlider3
f = figure;
ax = axes(f);
ax.Units = 'pixels';
ax.Position = [75 75 325 280]
c = uicontrol;
c.String = 'Plot Data';
c.Callback = @plotButtonPushed;

    function plotButtonPushed(src,event)
        bar(randn(1,5));
    end

end