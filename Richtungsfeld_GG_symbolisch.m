%Parameter, die in den DGLs variiert werden können
o = 0.3;
r1 = 0.6;
r2 = 0.5;
mu = 0.1

%Definition der DGL-Systeme in symbolischer Weise, hier in Vektorform. Ein
%Vektoreintrag entspricht einer Gleichung/Zeile des Systems
syms t dgl1(x,y)
dgl1(x,y) = [y; y*(-o)+x*(-r1)];
dgl2(x,y) = [y; mu*(1-x^2)*y-x];

%findet Gleichgewichtspunkte durch Nullsetzen der DGL
S1 = solve(dgl1(x,y)==0, [x y], 'ReturnConditions', true);
S2 = solve(dgl2(x,y)==0, [x y], 'ReturnConditions', true);

%gibt die Koordinaten der Gleichgewichtspunkte (bis zu fünf), setzt dabei
%die Werte -2 bis 2 für eventuelle Parameter ein
solx1 = subs(S1.x, S1.parameters, [-2 -1 0 1 2]);
soly1 = subs(S1.y, S1.parameters, [-2 -1 0 1 2]);

solx2 = subs(S2.x, S2.parameters, [-2 -1 0 1 2]);
soly2 = subs(S2.y, S2.parameters, [-2 -1 0 1 2]);

%erstellt Gitter, setzt Punkte in DGLs ein
[X,Y] = meshgrid(-10:1:10,-10:1:10);
d1=dgl1(X,Y);
d2=dgl2(X,Y);

%normieren (durch die geschweiften Klammern wird die erste bzw. zweite
%Gleichung des Systems aufgerufen)
dX1n = d1{1}./sqrt(d{1}.^2+d1{2}.^2);
dY1n = d1{2}./sqrt(d1{1}.^2+d1{2}.^2);

dX2n = d2{1}./sqrt(d2{1}.^2+d2{2}.^2);
dY2n = d2{2}./sqrt(d2{1}.^2+d2{2}.^2);

%wandelt symbolische Funktion in functionhandle um, weil ode45 nur damit
%arbeitet
fun1 = matlabFunction(dgl1,'Vars',{t,[x;y]});
fun2 = matlabFunction(dgl2,'Vars',{t,[x;y]});

%löst DGL, zweites Argument sind Werte für t, drittes die Startkoordinaten.
%Durch Ausprobieren muss bestimmt werden, welche Startkoordinaten ein gutes
%Bild ergeben
[t1,y1] = ode45(fun1,[0 20],[0 9]);
[t2,y2] = ode45(fun2,[0 20],[2 0]);


figure;
%erste DGL
%Richtungsfeld
subplot(1,2,1),q=quiver(X,Y,dX1n,dY1n);
q.Color = '#DEDEDE';
q.ShowArrowHead = 'off';
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx1, soly1, 'o');
hold on
%Lösungskurve
plot(y1(:,1),y1(:,2))

%zweite DGL
%Richtungsfeld
subplot(1,2,2),q2=quiver(X,Y,dX2n,dY2n);
q2.Color = '#DEDEDE';
q2.ShowArrowHead = 'off';
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx2, soly2, 'o');
hold on
%Lösungskurve
plot(y2(:,1),y2(:,2))





