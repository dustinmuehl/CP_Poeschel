%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : main_code.m                                                   %
%                                                                         %
% Authors : Moritz Amann , Dustin Mühlhäuser, Ilya Shapiro                % 
%           Benedikt Leibinger, Isabell Giers                             %   
% Date    : 22.12.2019                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Aufräumen
clc
clear
close all

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



%symbolische Bildung der zugehörigen Jacobi-Matrizen
A = jacobian(dgl1, [x, y]);
B = jacobian(dgl2, [x, y]);
C = jacobian(dgl3, [x, y]);
D = jacobian(dgl4, [x, y]);

%zugehörige EV und EW der Jacobi für Klassifikation
[VA, DA] = eig(A);
[VB, DB] = eig(B);
[VC, DC] = eig(C);
[VD, DD] = eig(D);

%Umwandlung DGL in functionhandle damit ode45 die Eingabe der DGL nutzen kann
syms t
fun1 = matlabFunction(dgl1,'Vars',{t,[x;y]});
fun2 = matlabFunction(dgl2,'Vars',{t,[x;y]});
fun3 = matlabFunction(dgl3,'Vars',{t,[x;y]});
fun4 = matlabFunction(dgl4,'Vars',{t,[x;y]});




%%
%Gleichgewichtspunkte

syms x y
dgl1(x,y)= dgl1;
dgl2(x,y)= dgl2;
dgl3(x,y)= dgl3;
dgl4(x,y)= dgl4;

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
for i = 1:size(solx1,1)
    Matrix = subs(A, x, solx1(i));
    [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx1(i), soly1(i));
    %Startwerte plotten
    if Farbe == 1       %konvergiert nach innen
        FC = 'r';
    elseif Farbe == -1  %divergiert nach ausssen
        FC = 'g';
    end
    
    for j = 1:size(AWerte,2) 
        [t,y] = ode45(fun1,[0 20*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
        plot(y(:,1),y(:,2), FC)
    end
    
    
    
    %Eigenraeume plotten falls notwendig
    if ERaum1 ~= 0
        plot(ERaum1(1,:),ERaum1(2,:), 'b','LineWidth',1.3)
    end
    if ERaum2 ~= 0
        plot(ERaum2(1,:),ERaum2(2,:), 'b','LineWidth',1.3)
    end
end




% %zweiter Plot
% subplot(2,2,2),q1=quiver(X,Y,dX2n,dY2n, 0.5);
% q1.Color = '#DEDEDE';
% q1.ShowArrowHead = 'off';
% title('DGL 2')
% axis([-10 10 -10 10])
% hold on
% %GG-Punkte
% plot(solx2, soly2, 'o');
% hold on
% %Loesungen und Eigenraume 
% %(fehlt noch: Spezifikation im Fall Sattel)
% for i = 1:size(solx2,1)
%     Matrix = subs(B, x, solx2(i));
%     [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx2(i), soly2(i));
%     %Startwerte plotten
%     if Farbe == 1       %konvergiert nach innen
%         FC = 'r';
%     elseif Farbe == -1  %divergiert nach ausssen
%         FC = 'g';
%     end
%     
%     for j = 1:size(AWerte,2) 
%         [t,y] = ode45(fun1,[0 20*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
%         plot(y(:,1),y(:,2), FC)
%     end
%     
%     
%     %Eigenraeume plotten falls notwendig
%     if ERaum1 ~= 0
%         plot(ERaum1(1,:),ERaum1(2,:))
%     end
%     if ERaum2 ~= 0
%         plot(ERaum2(1,:),ERaum2(2,:))
%     end
% end



% %dritter Plot
% subplot(2,2,3),q1=quiver(X,Y,dX1n,dY1n, 0.5);
% q1.Color = '#DEDEDE';
% q1.ShowArrowHead = 'off';
% title('DGL 1')
% axis([-10 10 -10 10])
% hold on
% %GG-Punkte
% plot(solx3, soly3, 'o');
% hold on
% %Loesungen und Eigenraume 
% %(fehlt noch: Spezifikation im Fall Sattel)
% for i = 1:size(solx3,1)
%     Matrix = subs(C, x, solx3(i));
%     [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx3(i), solx3(i));
%     %Startwerte plotten
%     if Farbe == 1       %konvergiert nach innen
%         FC = 'r';
%     elseif Farbe == -1  %divergiert nach ausssen
%         FC = 'g';
%     end
%     
%     for j = 1:size(AWerte,2) 
%         [t,y] = ode45(fun1,[0 20*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
%         plot(y(:,1),y(:,2), FC)
%     end
%     
%     
%     %Eigenraeume plotten falls notwendig
%     if ERaum1 ~= 0
%         plot(ERaum1(1,:),ERaum1(2,:))
%     end
%     if ERaum2 ~= 0
%         plot(ERaum2(1,:),ERaum2(2,:))
%     end
% end



% %vierter Plot
% subplot(2,2,4),q1=quiver(X,Y,dX4n,dY4n, 0.5);
% q1.Color = '#DEDEDE';
% q1.ShowArrowHead = 'off';
% title('DGL 4')
% axis([-10 10 -10 10])
% hold on
% %GG-Punkte
% plot(solx4, soly4, 'o');
% hold on
% %Loesungen und Eigenraume 
% %(fehlt noch: Spezifikation im Fall Sattel)
% for i = 1:size(solx4,1)
%     Matrix = subs(D, x, solx4(i));
%     [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx4(i), soly4(i));
%     %Startwerte plotten
%     [t,y] = ode45(fun4,[0 20*Farbe],Anfangswerte);
%     
%     if Farbe == 1       %konvergiert nach innen
%         FC = 'r';
%     elseif Farbe == -1  %divergiert nach ausssen
%         FC = 'g';
%     end
%     
%     plot(y(:,1),y(:,2), FC)
%     
%     %Eigenraeume plotten falls notwendig
%     if ERaum1 ~= 0
%         plot(ERaum1(1,:),ERaum1(2,:))
%     end
%     if ERaum2 ~= 0
%         plot(ERaum2(1,:),ERaum2(2,:))
%     end
% end



%aus Uebersichtsgruenden erstmal weggelassen
%Code ist identisch, einfach kopieren nur quiver anpassen, A und solx1,soly1 austauschen



%%
%Funktionen

%function [Anfangswerte, ERaum1,ERaum2]=fct_Klassifikation(J,xvalue,yvalue)


%entartete Fälle, dh doppelte EW




%entartete Knoten, dh det=0



%reele EW
function [A,B1,B2,R] = fkt_Sattel(VJ,WJ,xValue,yValue,r,k)
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


% function [A,B1,B2]=fkt_Knoten(VJ,WJ,xValue,yValue,k,r)
% end

function [A,B1,B2] = fkt_Knoten(J,xValue,yValue,k,r)
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
function [A] = fkt_Zentrum(xValue,yValue,a)
    a1 = [xValue+a yValue];
    a2 = [xValue+2*a yValue];
    a3 = [xValue+3*a yValue];
    a4 = [xValue+4*a yValue];
    a5 = [xValue+5*a yValue];
    A = [a1; a2; a3; a4; a5];
end

function [A] = fkt_Strudel(xValue,yValue,r)
        A = [xValue+r yValue];
end

function [A, B1, B2] = fkt_lin_Knoten(xWert, yWert, EV0)
    B1 = [xWert+3*EV0(1) xWert-3*EV0(1); yWert+3*EV0(2) yWert-3*EV0(2)];
    B2 = 0; %nur ein Eigenraum
    
    if EV0(1) == 0
        k = 0.1;
    else
        k = 0;
    end
        
    a1 = [xWert-2*EV0(1)+k; yWert-2*EV0(2)+0.1];
    a2 = [xWert-2*EV0(1)-k; yWert-2*EV0(2)-0.1];
    a3 = [xWert-1*EV0(1)+k; yWert-1*EV0(2)+0.1];
    a4 = [xWert-1*EV0(1)-k; yWert-1*EV0(2)-0.1];
    a5 = [xWert+1*EV0(1)+k; yWert+1*EV0(2)+0.1];
    a6 = [xWert+1*EV0(1)-k; yWert+2*EV0(2)-0.1];
    a7 = [xWert+2*EV0(1)+k; yWert+2*EV0(2)+0.1];
    a8 = [xWert+2*EV0(1)-k; yWert+2*EV0(2)-0.1];
    
    A = [a1 a2 a3 a4 a5 a6 a7 a8];
end


function [A, B1, B2] = fkt_Stern(xWert,yWert)
    a1 = [xWert+0.1; yWert+0];
    a2 = [xWert-0.1; yWert+0];
    a3 = [xWert+0; yWert+0.1];
    a4 = [xWert+0; yWert-0.1];
    a5 = [xWert+0.1; yWert+0.1];
    a6 = [xWert-0.1; yWert+0.1];
    a7 = [xWert+0.1; yWert-0.1];
    a8 = [xWert-0.1; yWert-0.1];

    A = [a1 a2 a3 a4 a5 a6 a7 a8];
    % kein Eigenraum
    B1 = 0;
    B2 = 0;
end

%%
%Klassifikation

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
                    %entartet, instab Knoten, nicht diagonalisierbar
                    [A, B1, B2] = fkt_Stern(xWert, yWert);
                    R = -1;
                else
                    %entartet, stab Knoten, diagonalisierbar
                    [A, B1, B2, R] = fkt_ent_Knoten(INPUTS);
                end
            else
                if rank(VJ) == 1
                    %entartet, stab Knoten, nicht diagonalisierbar
                    [A, B1, B2] = fkt_Stern(xWert, yWert);
                    R = 1;
                else
                    %entartet, stab Knoten, diagonalisierbar
                    [A, B1, B2] = fkt_ent_Knoten(INPUTS);
                    R = 1;
                end
            end
        end
    
       
    
    else %lambda ungleich mu
        if lambda * mu == 0
            if mu == 0
                if lambda > 0
                    %lineare Expansion mit mu = 0
                    [A, B1, B2] = fkt_lin_Knoten(xWert, yWert, [VJ(1,2); VJ(2,2)]);
                    R = -1;
                else
                    %lineare Kontraktion
                    [A, B1, B2] = fkt_lin_Knoten(xWert, yWert, [VJ(1,2); VJ(2,2)]);
                    R = 1;
                end
            else
                if mu > 0
                    %lineare Expansion mit lambda=0
                    [A, B1, B2] = fkt_lin_Konten(xWert, yWert, [VJ(1,1); VJ(2,1)]);
                    R = -1;
                else 
                    %lineare Kontraktion
                    [A, B1, B2] = fkt_lin_Knoten(xWert, yWert, [VJ(1,1); VJ(2,1)]);
                    R = 1;
                end
            end
            
        else
            if abs(imag(lambda)) < 0.001
                if lambda * mu < 0
                    [A, B1, B2, R] = fkt_Sattel(VJ,DJ,xWert,yWert,2,5);
                else
                    if lambda > 0
                        %instabil
                        [A, B1, B2] = fkt_Knoten(VJ,DJ,xWert,yWert,2,5);
                        R = -1;
                    else
                        %stabil
                        [A, B1, B2] = fkt_Knoten(VJ,DJ,xWert,yWert,2,5);
                        R = 1;
                    end
                end
                
            else %lambda imaginaer
                if abs(real(lambda)) < 0.001
                    A = fkt_Zentrum(xWert,yWert,1);
                    R = 1;
                    B1 = 0;
                    B2 = 0;
                else
                    if real(lambda) > 0
                        %instabil
                        A = fkt_Strudel(xWert,yWert,0.05);
                        R = -1;
                        B1 = 0;
                        B2 = 0;
                    else
                        %stabil
                        A = fkt_Strudel(xWert,yWert,0.05);
                        R = 1;
                        B1 = 0;
                        B2 = 0;
                    end
                end
            end
        end
    end

end



%%
%Video