%Parameter, die in den DGLs variiert werden können
p1 = 2;
w1 = 9;
p2 = 1;

% Jacobi-Matrizen der Systeme
A = [0 1; -w1 -p1];
B1 = [0 1, -1 -p2];
B2 = [0 1, +1 -p2];
% C = [];
% E = [];

%VA ist Matrix der Eigenvektoren zu A, DA ist Matrix der Eigenwerte zu A
[VA,DA] = eig(A);
% [VB1,DB1] = eig(B1)
% [VB2,DB2] = eig(B1)
% [VC, DC] = eig(C)
% [VE, BE] = eig(E)


%Definition der DGL-Systeme in symbolischer Weise, hier in Vektorform. Ein
%Vektoreintrag entspricht einer Gleichung/Zeile des Systems
syms t dgl1(x,y)
dgl1(x,y) = [y; y*(-p1)+x*(-w1)];
dgl2(x,y) = [y; -sin(x) - p2 *y];

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
dX1n = d1{1}./sqrt(d1{1}.^2+d1{2}.^2);
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
[t1,y1] = ode45(fun1,[0 20],[-4 6]);
[t2,y2] = ode45(fun2,[0 20],[2 0]);


figure;
%erste DGL
%Richtungsfeld
subplot(2,2,1),q=quiver(X,Y,dX1n,dY1n, 0.5);
q.Color = '#DEDEDE';
q.ShowArrowHead = 'off';

axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx1, soly1, 'o');
hold on
%Lösungskurve
plot(y1(:,1),y1(:,2))
% �berpr�fen ob Eigenwerte reell sind
if abs(imag(DA(1,1))) < 0.001
    % Plotten der Eigenvektoren an Gleichgewichtspunkt
    plot([solx1-3*VA(1,1) solx1+ 3*VA(1,1)],[soly1-3*VA(2,1) soly1+3*VA(2,1)])
end
if abs(imag(DA(2,2))) < 0.001
    plot([solx1-3*VA(1,2) solx1+ 3*VA(1,2)],[soly1-3*VA(2,2) soly1+3*VA(2,2)])
end

%zweite DGL
%Richtungsfeld
subplot(2,2,2),q2=quiver(X,Y,dX2n,dY2n, 0.5);
q2.Color = '#DEDEDE';
q2.ShowArrowHead = 'off';
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx2, soly2, 'o');
hold on
%Lösungskurve
plot(y2(:,1),y2(:,2))



% Platzhalter f�r andere Systeme

subplot(2,2,3),q2=quiver(X,Y,dX2n,dY2n, 0.5);
q2.Color = '#DEDEDE';
q2.ShowArrowHead = 'off';
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx2, soly2, 'o');
hold on
%Lösungskurve
plot(y2(:,1),y2(:,2))

subplot(2,2,4),q2=quiver(X,Y,dX2n,dY2n, 0.5);
q2.Color = '#DEDEDE';
q2.ShowArrowHead = 'off';
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx2, soly2, 'o');
hold on
%Lösungskurve
plot(y2(:,1),y2(:,2))




