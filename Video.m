%Input range festlegen 
z = 0:0.05:10; 
y = sin(2*z);
x = cos(2*z);

%Erzeugung des Videos 
v = VideoWriter('DGL.avi','Motion JPEG AVI');

%Oeffnen des Videos 
open(v);

%Animation des Videos und Aufnahme 
curve = animatedline('LineWidth', 2);
set(gca, 'XLim', [-1.5 1.5], 'YLim', [-1.5 1.5], 'ZLim', [0 10]); 
view(43, 24);
hold on;

for i = 1:length(z)
    addpoints(curve, x(i), y(i), z(i));
    drawnow
    pause(0.01);
    currFrame = getframe;
    writeVideo(v,currFrame);
end 

%Beenden des Videos 
close(v);