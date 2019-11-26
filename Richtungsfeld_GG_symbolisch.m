

o = 0.5;
r1 = 0.3;
r2 = 0.4;

syms P(x,y) Q(x,y) R(x,y)
%DGLs symbolisch angeben
P(x,y) = y;                         %erste Zeile math. Pendel linearisiert
Q(x,y) = y*(-o)+x*(-r1);
R(x,y) = -sin(x)+y*(-r2);

%ermittelt Gleichgewichtspunkte
S1 = solve([P(x,y)==0, Q(x,y)==0], [x,y], 'ReturnConditions', true);
S2 = solve([P(x,y)==0, R(x,y)==0], [x,y], 'ReturnConditions', true);

%gibt x-, y-Koordinaten von bis zu fünf Gleichgewichtspunkten
S2.x
solx1 = subs(S1.x, S1.parameters, [-2 -1 0 1 2]);
soly1 = subs(S1.y, S1.parameters, [-2 -1 0 1 2]);

solx2 = subs(S2.x, S2.parameters, [-2 -1 0 1 2]);
soly2 = subs(S2.y, S2.parameters, [-2 -1 0 1 2]);

%für Richtungsfeld
[X,Y] = meshgrid(-10:1:10,-10:1:10);
dX = P(X,Y);
dY1 = Q(X,Y);
dY2 = R(X,Y);

%normieren
dY1n = dY1./sqrt(dX.^2+dY1.^2);
dX1n = dX./sqrt(dX.^2+dY1.^2);

dY2n = dY2./sqrt(dX.^2+dY2.^2);
dX2n = dX./sqrt(dX.^2+dY2.^2);


figure;
%Richtungsfeld
subplot(1,2,1),quiver(X,Y,dX1n,dY1n)
axis tight
hold on
%GG-Punkte
plot(solx1, soly1, 'o');
subplot(1,2,2),quiver(X,Y,dX2n,dY2n)
axis tight
hold on
plot(solx2, soly2, 'o');







