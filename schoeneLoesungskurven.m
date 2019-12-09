%Parameter, die in den DGLs variiert werden koennen
%z.B.: Sattel p1=2,w1=1; stabiler Strudel 0.5,1; instabiler Strudel -0.5,1
p1 = 2;
w1 = 1;
p2 = 0.5;
m1=2;
m2=1;

% Jacobi-Matrizen der Systeme
A = [0 1; -w1 -p1];
B1 = [0 1; -1 -p2];
B2 = [0 1; +1 -p2];
syms z
C = [m1-3*z^2 0; 0 -1];
E = [m2^3 - 3*m2^2*z*(2+z) + 3*m2*(-1+z^2) + z*18 + 9*z - 12*z^2 - 5*z^3 0; 0 -1];

%VA ist Matrix der Eigenvektoren zu A, DA ist Matrix der Eigenwerte zu A
[VA,DA] = eig(A);
[VB1,DB1] = eig(B1);
[VB2,DB2] = eig(B1);
% [VC, DC] = eig(C)
% [VE, BE] = eig(E)


%Definition der DGL-Systeme in symbolischer Weise, hier in Vektorform. Ein
%Vektoreintrag entspricht einer Gleichung/Zeile des Systems
syms t dgl1(x,y)
dgl1(x,y) = [y; y*(-p1)+x*(-w1)];
dgl2(x,y) = [y; -sin(x) - p2 *y];
dgl3(x,y) = [y; m1*x-x^3];
dgl4(x,y) = [y; -x*(x^2+(m2)^2-3)*(x^2+3*x-m2)];

%findet Gleichgewichtspunkte durch Nullsetzen der DGL
S1 = solve(dgl1(x,y)==0, [x y], 'ReturnConditions', true);
S2 = solve(dgl2(x,y)==0, [x y], 'ReturnConditions', true);
S3 = solve(dgl3(x,y)==0, [x y], 'ReturnConditions', true);
S4 = solve(dgl4(x,y)==0, [x y], 'ReturnConditions', true);

%gibt die Koordinaten der Gleichgewichtspunkte (bis zu fuenf), setzt dabei
%die Werte -2 bis 2 fuer eventuelle Parameter ein
solx1 = subs(S1.x, S1.parameters, [-2 -1 0 1 2]);
soly1 = subs(S1.y, S1.parameters, [-2 -1 0 1 2]);

solx2 = subs(S2.x, S2.parameters, [-2 -1 0 1 2]);
soly2 = subs(S2.y, S2.parameters, [-2 -1 0 1 2]);

solx3 = subs(S3.x, S3.parameters, [-2 -1 0 1 2]);
soly3 = subs(S3.y, S3.parameters, [-2 -1 0 1 2]);

solx4 = subs(S4.x, S4.parameters, [-2 -1 0 1 2]);
soly4 = subs(S4.y, S4.parameters, [-2 -1 0 1 2]);

%erstellt Gitter, setzt Punkte in DGLs ein
[X,Y] = meshgrid(-10:1:10,-10:1:10);
d1=dgl1(X,Y);
d2=dgl2(X,Y);
d3=dgl3(X,Y);
d4=dgl4(X,Y);

%normieren (durch die geschweiften Klammern wird die erste bzw. zweite
%Gleichung des Systems aufgerufen)
dX1n = d1{1}./sqrt(d1{1}.^2+d1{2}.^2);
dY1n = d1{2}./sqrt(d1{1}.^2+d1{2}.^2);

dX2n = d2{1}./sqrt(d2{1}.^2+d2{2}.^2);
dY2n = d2{2}./sqrt(d2{1}.^2+d2{2}.^2);

dX3n = d3{1}./sqrt(d3{1}.^2+d3{2}.^2);
dY3n = d3{2}./sqrt(d3{1}.^2+d3{2}.^2);

dX4n = d4{1}./sqrt(d4{1}.^2+d4{2}.^2);
dY4n = d4{2}./sqrt(d4{1}.^2+d4{2}.^2);

%wandelt symbolische Funktion in functionhandle um, weil ode45 nur damit
%arbeitet
fun1 = matlabFunction(dgl1,'Vars',{t,[x;y]});
fun2 = matlabFunction(dgl2,'Vars',{t,[x;y]});
fun3 = matlabFunction(dgl3,'Vars',{t,[x;y]});
fun4 = matlabFunction(dgl4,'Vars',{t,[x;y]});


figure;
%erste DGL
%Richtungsfeld
subplot(2,2,1),q1=quiver(X,Y,dX1n,dY1n, 0.5);
q1.Color = '#DEDEDE';
q1.ShowArrowHead = 'off';
title('DGL 1')
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx1, soly1, 'o');
hold on

% Ueberpruefen ob Eigenwerte reell sind
sz1 = size(solx1);
if abs(imag(DA(1,1))) < 0.001
    % Plotten der Eigenvektoren an Gleichgewichtspunkt (Fall Sattel)
    for i = 1:sz1(2)
        plot([solx1(i)-2*VA(1,1) solx1(i)+ 2*VA(1,1)],[soly1(i)-2*VA(2,1) soly1(i)+2*VA(2,1)])
    end
    %errechnen und plotten der Loesungskurven
    [t1,y1] = ode45(fun1,[0 20],[-5 10]);
    plot(y1(:,1),y1(:,2))
    [t1,y1] = ode45(fun1,[0 20],[5 10]);
    plot(y1(:,1),y1(:,2))
    [t1,y1] = ode45(fun1,[0 20],[-5 -10]);
    plot(y1(:,1),y1(:,2))
    [t1,y1] = ode45(fun1,[0 20],[5 -10]);
    plot(y1(:,1),y1(:,2))

elseif real(DA(1,1)) > 0
    %Fall instabiler Strudel
    [t1,y1] = ode45(fun1,[0 20],[double(solx1(1))+0.1 double(soly1(1))]);
    plot(y1(:,1),y1(:,2))

elseif real(DA(1,1)) < 0
    %Fall stabiler Strudel
    [t1,y1] = ode45(fun1,[0 20],[double(solx1(1)) double(soly1(1))+10]);
    plot(y1(:,1),y1(:,2))
end

if abs(imag(DA(2,2))) < 0.001
    for i = 1:sz1(2)
        plot([solx1(i)-2*VA(1,2) solx1(i)+ 2*VA(1,2)],[soly1(i)-2*VA(2,2) soly1(i)+2*VA(2,2)])
    end
end

%zweite DGL
%Richtungsfeld
subplot(2,2,2),q2=quiver(X,Y,dX2n,dY2n, 0.5);
q2.Color = '#DEDEDE';
q2.ShowArrowHead = 'off';
title('DGL 2')
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx2, soly2, 'o');
hold on

sz2 = size(solx2);
if abs(imag(DB1(1,1))) < 0.001
    % Plotten der Eigenvektoren an Gleichgewichtspunkt
    for i = 1:sz2(2)
        plot([solx2(i)-2*VB1(1,1) solx2(i)+ 2*VB1(1,1)],[soly2(i)-2*VB1(2,1) soly2(i)+2*VB1(2,1)])
    end
    %Loesungskurve
    [t2,y2] = ode45(fun2,[0 20],[-3 3]);
    plot(y2(:,1),y2(:,2))
else
    %Fall imaginaere Eigenwerte
    [tim,yim] = ode45(fun2,[0 20],[double(solx2(2)) double(soly2(2))+0.0001]);
    plot(yim(:,1),yim(:,2))
    [tim,yim] = ode45(fun2,[0 20],[double(solx2(2)) double(soly2(2))-0.0001]);
    plot(yim(:,1),yim(:,2))
    [tim,yim] = ode45(fun2,[0 20],[double(solx2(4)) double(soly2(4))+0.0001]);
    plot(yim(:,1),yim(:,2))
    [tim,yim] = ode45(fun2,[0 20],[double(solx2(4)) double(soly2(4))-0.0001]);
    plot(yim(:,1),yim(:,2))
end
if abs(imag(DB1(2,2))) < 0.001
    for i = 1:sz2(2)
        plot([solx2(i)-2*VB1(1,2) solx2(i)+ 2*VB1(1,2)],[soly2(i)-2*VB1(2,2) soly2(i)+2*VB1(2,2)])
    end
end





%dritte DGL

subplot(2,2,3),q3=quiver(X,Y,dX3n,dY3n, 0.5);
q3.Color = '#DEDEDE';
q3.ShowArrowHead = 'off';
axis([-10 10 -10 10])
title('DGL3')
hold on
%GG-Punkte
plot(solx3, soly3, 'o');
hold on
%Loesungskurve
[t3,y3] = ode45(fun3,[0 20],[2 0]);
plot(y3(:,1),y3(:,2))


%vierte DGL

subplot(2,2,4),q4=quiver(X,Y,dX4n,dY4n, 0.5);
q4.Color = '#DEDEDE';
q4.ShowArrowHead = 'off';
title('DGL4')
axis([-10 10 -10 10])
hold on
%GG-Punkte
plot(solx4, soly4, 'o');
hold on
%Loesungskurve
[t4,y4] = ode45(fun4,[0 20],[2 0]);
plot(y4(:,1),y4(:,2))
