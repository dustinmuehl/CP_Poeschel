%Lege die Konstanten Fest
o = 0.3;
r1 = 0.6;
r2 = 0.5;

%erzeuge Gitter
[y1,y2] = meshgrid(-10:1:10,-10:1:10);
%Berechne Werte-Arrays der auf dgl System 1. Grades umgeformten dgl 
dy1 = y2;
dy2 = y1.*(-o) + y2.*(-r1) ;
%normiere quiver Pfeile
dy2u = dy2./sqrt(dy1.^2+dy2.^2);
dy1u = dy1./sqrt(dy1.^2+dy2.^2);

dx1 = y2;
dx2 = - sin(y1) + y2.*(-r2) ;
dx2u = dx2./sqrt(dx1.^2+dx2.^2);
dx1u = dx1./sqrt(dx1.^2+dx2.^2);

%Erzeuge Plotfenster und verwende quiver-Funktion um Richtungsfeld zu
%zeichnen
figure;
subplot(1,2,1),quiver(y1,y2,dy1u,dy2u);
subplot(1,2,2),quiver(y1,y2,dx1u,dx2u);