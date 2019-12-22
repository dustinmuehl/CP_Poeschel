%%
%Input
%Input Parameter
p1 = 2;
w1 = 1;
p2 = 3;
m1 = 2;
m2 = 1;

%Input DGL´s
syms x y
dgl1 = [y; y*(-p1)+x*(-w1)];
dgl2 = [y; -sin(x) - p2 *y];
dgl3 = [m1*x-x^3; -y];
dgl4 = [-x*(x^2+(m2)^2-3)*(x^2+3*x-m2); -y];

%%

%symbolische Bildung der zugehörigen Jacobi-Matrizen
A = jacobian(dgl1, [x, y]);
B = jacobian(dgl2, [x, y]);
C = jacobian(dgl3, [x, y])
D = jacobian(dgl4, [x, y]);

%zugehörige EV und EW der Jacobi für Klassifikation
[VA, DA] = eig(A);
[VB, DB] = eig(B);
[VC, DC] = eig(C);
[VD, DD] = eig(D);

%Umwandlung DGL in functionhandle damit ode45 die Eingabe der DGL nutzen kann
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

function [A, B1, B2, R] = fkt_Klassifikation(J, xWert, yWert)

    [VJ,DJ] = eig(J);
    lambda = DJ(1,1);
    mu = DJ(2,2);
    
    if mu == lambda
        if mu * lambda == 0
            [A, B1, B2, R] = fkt_entar_Knoten_EW0(INPUTS);
        else
            if lambda > 0
                if rank(VJ) == 1
                    [A, B1, B2, R] = fkt_entar_instab_Knoten_Ndiag(INPUTS);
                else
                    [A, B1, B2, R] = fkt_entar_instab_Knoten_diag(INPUTS);
                end
            else
                if rank(VJ) == 1
                    [A, B1, B2, R] = fkt_entar_stab_Knoten_Ndiag(INPUTS);
                else
                    [A, B1, B2, R] = fkt_entar_stab_Knoten_diag(INPUTS);
                end
            end
        end
    
        
        
    else %lambda ungleich mu
        if lambda * mu == 0
            if mu == 0
                if lambda > 0
                    [A, B1, B2, R] = fkt_lin_Expansion_mu0(INPUTS);
                else
                    [A, B1, B2, R] = fkt_lin_Kontraktion_mu0(INPUTS);
                end
            else
                if mu > 0
                    [A, B1, B2, R] = fkt_lin_Expansion_lambda0(INPUTS);
                else 
                    [A, B1, B2, R] = fkt_lin_Kontraktion_lambda0(INPUTS);
                end
            end
            
        else
            if abs(imag(lambda)) < 0.001
                if lambda * mu < 0
                    [A, B1, B2, R] = fkt_Sattel(INPUTS);
                else
                    if lambda > 0
                        [A, B1, B2] = fkt_Knoten(INPUTS);
                        R = -1; %instabil
                    else
                        [A, B1, B2] = fkt_Knoten(INPUTS);
                        R = 1; %stabil
                    end
                end
                
            else %lambda imaginaer
                if abs(real(lambda)) < 0.001
                    A = fkt_Zentrum(INPUTS);
                    R = 1;
                    B1 = 0;
                    B2 = 0;
                else
                    if real(lambda) > 0
                        A = fkt_Strudel(INPUTS);
                        R = -1; %instabil
                        B1 = 0;
                        B2 = 0;
                    else
                        A = fkt_Strudel(INPUTS);
                        R = 1; %stabil
                        B1 = 0;
                        B2 = 0;
                    end
                end
            end
        end
    end

end


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


%entartete Fälle, dh doppelte EW




%entartete Knoten, dh det=0



%reele EW
function [A,B1,B2,R] = fkt_sattel(VJ,WJ,xValue,yValue,r,k)
    a1 = [xValue+k*VJ(1,1)+r yValue+k*VJ(2,1)+r];
    a2 = [xValue+k*VJ(1,1)+r yValue+k*VJ(2,1)-r];
    a3 = [xValue+k*VJ(1,1)-r yValue+k*VJ(2,1)+r];
    a4 = [xValue+k*VJ(1,1)-r yValue+k*VJ(2,1)-r];
    a5 = [xValue-k*VJ(1,1)+r yValue-k*VJ(2,1)+r];
    a6 = [xValue-k*VJ(1,1)+r yValue-k*VJ(2,1)-r];
    a7 = [xValue-k*VJ(1,1)-r yValue-k*VJ(2,1)+r];
    a8 = [xValue-k*VJ(1,1)-r yValue-k*VJ(2,1)-r];
    A = [a1; a2; a3; a4; a5; a6; a7; a8];
    B1 = VJ(:,1); 
    B2 = VJ(:,2);
if WJ(1,1) < WJ(2,2)
    R=1;
else
    R=-1;
end   
end

function [A,B1,B2]=fkt_knoten(VJ,WJ,xValue,yValue,k,r)
end

function [A,B1,B2] = fct_Knoten(J,xValue,yValue,k,r)
% J ist die Jacobi Determinante
% xValue und yValue sind die Koordinaten des GG-Punkts
% k ist der Wert, mit dem der Eigenvektor multipliziert wird (guter Wert: 5)
% r ist der Wert, um den der Punkt dann verschoben wird (guter Wert: 2)
%
% Ausgabe: 8x2 Matrix, wobei eine Zeile aus x und y-Wert des Startpunkts
%           besteht

    [VJ,WJ] = eig(J);
        if abs(imag(WJ(1,1))) < 0.001
            if abs(WJ(1,1)) < abs(WJ(2,2))
                a1 = [xValue+k*VJ(1,1)+r yValue+k*VJ(1,2)+r];
                a2 = [xValue+k*VJ(1,1)+r yValue+k*VJ(1,2)-r];
                a3 = [xValue+k*VJ(1,1)-r yValue+k*VJ(1,2)+r];
                a4 = [xValue+k*VJ(1,1)-r yValue+k*VJ(1,2)-r];
                a5 = [xValue-k*VJ(1,1)+r yValue-k*VJ(1,2)+r];
                a6 = [xValue-k*VJ(1,1)+r yValue-k*VJ(1,2)-r];
                a7 = [xValue-k*VJ(1,1)-r yValue-k*VJ(1,2)+r];
                a8 = [xValue-k*VJ(1,1)-r yValue-k*VJ(1,2)-r];
            else
                a1 = [xValue+k*VJ(2,1)+r yValue+k*VJ(2,2)+r];
                a2 = [xValue+k*VJ(2,1)+r yValue+k*VJ(2,2)-r];
                a3 = [xValue+k*VJ(2,1)-r yValue+k*VJ(2,2)+r];
                a4 = [xValue+k*VJ(2,1)-r yValue+k*VJ(2,2)-r];
                a5 = [xValue-k*VJ(2,1)+r yValue-k*VJ(2,2)+r];
                a6 = [xValue-k*VJ(2,1)+r yValue-k*VJ(2,2)-r];
                a7 = [xValue-k*VJ(2,1)-r yValue-k*VJ(2,2)+r];
                a8 = [xValue-k*VJ(2,1)-r yValue-k*VJ(2,2)-r];
            end
            A = [a1; a2; a3; a4; a5; a6; a7; a8];
        end
    B1 = [xValue-3*VJ(1,1) xValue+3*VJ(1,1) ; yValue-3*VJ(1,2) yValue+3*VJ(1,2)];
    B2 = [xValue+3*VJ(2,1) xValue-3*VJ(2,1) ; yValue+3*VJ(2,2) yValue-3*VJ(2,2)];
end


% komplexe EW
function [A] = fkt_zentrum(xValue,yValue,a)
    a1 = [xValue+a yValue];
    a2 = [xValue+2*a yValue];
    a3 = [xValue+3*a yValue];
    a4 = [xValue+4*a yValue];
    a5 = [xValue+5*a yValue];
    A = [a1; a2; a3; a4; a5];
end

function [A] = fkt_strudel(xValue,yValue,r)
        A = [xValue+r yValue];
end



%%
%Video