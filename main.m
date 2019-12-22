%%
%Input 

%Input DGL´s
syms t
dgl1(x,y) = [y; y*(-p1)+x*(-w1)];
dgl2(x,y) = [y; -sin(x) - p2 *y];
dgl3(x,y) = [y; m1*x-x^3];
dgl4(x,y) = [y; -x*(x^2+(m2)^2-3)*(x^2+3*x-m2)];

%Input Jacobi-Matrizen
A=[]
B=[]
C=[]
D=[]

%wandelt symbolische Funktion in functionhandle um, weil ode45 nur damit
%arbeitet
fun1 = matlabFunction(dgl1,'Vars',{t,[x;y]});
fun2 = matlabFunction(dgl2,'Vars',{t,[x;y]});
fun3 = matlabFunction(dgl3,'Vars',{t,[x;y]});
fun4 = matlabFunction(dgl4,'Vars',{t,[x;y]});



%%
%Gleichgewichtspunkte

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



%%
%Klassifikation

[Anfangswerte,ERaum1,V2]=fct_Klassifikation(A,solx1,soly1);
[MB,VB1,VB2]=fct_Klassifikation(B,solx2,soly2);
[MC,VC1,VC2]=fct_Klassifikation(C,solx3,soly3);
[MD,VD1,VD2]=fct_Klassifikation(D,solx4,soly4);


%%
%Plots


%erster Plot
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
%Loesungen und Eigenraume 
%(fehlt noch: Spezifikation im Fall Sattel)
for i=1:size(solx1,2)
    [Anfangswerte,ERaum1,ERaum2]=fct_Klassifikation(A,solx1(i),soly1(i));
    
    if ERaum1==1 || ERaum1==-1 %Strudel
        [t,y] = ode45(fun1,[0 20*ERaum1],Anfangswerte);
        plot(y(:,1),y(:,2), 'g')
    elseif ERaum1==0 %Zentrum
        for j=1:size(Anfangswerte,1) %hier werden mehrere Lösungen eingezeichnet
            [t,y] = ode45(fun1,[0 20],Anfangswerte(j,:));
            plot(y(:,1),y(:,2), 'g')
        end
    else %Knoten (oder Sattel?)
        for j=1:size(Anfangswerte,1) %hier werden mehrere Lösungen eingezeichnet
            [t,y] = ode45(fun1,[0 20],Anfangswerte(j,:));
            plot(y(:,1),y(:,2), 'g')
        end
        %hier gibt es auch Eigenräume
        plot(ERaum1(1,:),ERaum1(2,:))
        plot(ERaum2(1,:),ERaum2(2,:))
    end
end

%zweiter Plot
subplot(2,2,2)
%dritter Plot
subplot(2,2,3)
%vierter Plot
subplot(2,2,4)

%aus Uebersichtsgruenden erstmal weggelassen
%Code ist identisch, einfach kopieren nur quiver anpassen, A und solx1,soly1 austauschen



%%
%Funktionen


%function [Anfangswerte, ERaum1,ERaum2]=fct_Klassifikation(J,xvalue,yvalue)



