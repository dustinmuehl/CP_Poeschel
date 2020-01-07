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
m1 = 3;
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
axis([-5 5 -5 5])
hold on
%GG-Punkte
%plot(solx1, soly1, 'o');
hold on

%Loesungen und Eigenraume 
for i = 1:size(solx1,1)
    Matrix = subs(A, x, solx1(i));
    [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx1(i), soly1(i));
    %Startwerte plotten
    if Farbe == -1       %konvergiert nach außen
        FC = 'r';
        PF = 'or';
    elseif Farbe == 1  %divergiert nach innen
        FC = 'g';
        PF = 'og';
    elseif Farbe == 1.0001 %neutral
        FC = 'b';
        PF = 'ob';
    end
    
    %Kreis um Loesung plotten
    plot(solx1(i), soly1(i), PF)
    
    if Farbe==0
    %wenn Jacobische 0, dann tue nichts    
    else
    %Eigenraeume plotten falls notwendig
    if ERaum1 ~= 0
        plot(ERaum1(1,:),ERaum1(2,:), 'k','LineWidth',1.0)
    end
    if ERaum2 ~= 0
        plot(ERaum2(1,:),ERaum2(2,:), 'k','LineWidth',1.0)
    end
    %Loesung plotten
    for j = 1:size(AWerte,2)
        [t,y] = ode45(fun1,[0 5*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
        plot(y(:,1),y(:,2), FC)
    end
    end
end

%zweiter Plot
subplot(2,2,2),q1=quiver(X,Y,dX2n,dY2n, 0.5);
q1.Color = '#DEDEDE';
q1.ShowArrowHead = 'off';
title('DGL 2')
axis([-10 10 -10 10])
hold on
%GG-Punkte
%plot(solx2, soly2, 'o');
hold on

%Loesungen und Eigenraume 
for i = 1:size(solx2,2)
    Matrix = subs(B, x, solx2(i));
    [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx2(i), soly2(i)); 
    
    %Startwerte plotten
    if Farbe == -1       %divergiert nach aussen
        FC = 'r';
        PF = 'or';
    elseif Farbe == 1  %konvergiert nach innen
        FC = 'g';
        PF = 'og';
    elseif Farbe == 1.0001 %neutral
        FC = 'b';
        PF = 'ob';
    end
    
    %Kreis um Loesung plotten
    plot(solx2(i), soly2(i), PF);

    if Farbe==0
    %wenn Jacobische 0, dann tue nichts    
    else
    %Eigenraeume plotten falls notwendig
    if ERaum1 ~= 0
        plot(ERaum1(1,:),ERaum1(2,:), 'k','LineWidth',1.0)
    end
    if ERaum2 ~= 0
        plot(ERaum2(1,:),ERaum2(2,:), 'k','LineWidth',1.0)
    end
    
    %Loesung plotten
    for j = 1:size(AWerte,2) 
        [t,y] = ode45(fun2,[0 5*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
        plot(y(:,1),y(:,2), FC)
    end
    end
end

%dritter Plot
subplot(2,2,3),q1=quiver(X,Y,dX3n,dY3n, 0.5);
q1.Color = '#DEDEDE';
q1.ShowArrowHead = 'off';
title('DGL 3')
axis([-5 5 -5 5])
hold on
%GG-Punkte
% plot(solx3, soly3, 'o');
hold on

%Loesungen und Eigenraume
for i = 1:size(solx3,1)
    Matrix = subs(C, x, solx3(i));
    [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx3(i), soly3(i));
    %Startwerte plotten
    if Farbe == -1       %divergiert nach aussen
        FC = 'r';
        PF = 'or';
    elseif Farbe == 1  %konvergiert nach innen
        FC = 'g';
        PF = 'og';
    elseif Farbe == 1.0001 %neutral
        FC = 'b';
        PF = 'ob';
    end
    
    %Kreis um Loesung plotten
    plot(solx3(i), soly3(i), PF)
    
    if Farbe==0
    %wenn Jacobische 0, dann tue nichts    
    else
    %Eigenraeume plotten falls notwendig
    if norm(ERaum1) ~= 0
        plot(ERaum1(1,:),ERaum1(2,:), 'k','LineWidth',1.0)
    end
    if norm(ERaum2) ~= 0
        plot(ERaum2(1,:),ERaum2(2,:), 'k','LineWidth',1.0)
    end
    
    %Loesung plotten
    for j = 1:size(AWerte,2)
        [t,y] = ode45(fun3,[0 5*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
        plot(y(:,1),y(:,2), FC)
    end
    end
end



%vierter Plot
subplot(2,2,4),q1=quiver(X,Y,dX4n,dY4n, 0.5);
q1.Color = '#DEDEDE';
q1.ShowArrowHead = 'off';
title('DGL 4')
axis([-5 5 -5 5])
hold on
%GG-Punkte
% plot(solx4, soly4, 'o');
hold on

%Loesungen und Eigenraume 
for i = 1:size(solx4,1)
    Matrix = subs(D, x, solx4(i));
    [AWerte, ERaum1, ERaum2, Farbe] = fkt_Klassifikation(Matrix, solx4(i), soly4(i));
    %Startwerte plotten
    if Farbe == -1       %divergiert nach aussen
        FC = 'r';
        PF = 'or';
    elseif Farbe == 1  %konvergiert nach innen
        FC = 'g';
        PF = 'og';
    elseif Farbe == 1.0001 %neutral
        FC = 'b';
        PF = 'ob';
    end
    
    %Kreis um Loesung plotten 
    plot(solx4(i), soly4(i), PF)

    if Farbe==0
    %wenn Jacobische 0, dann tue nichts    
    else
    %Eigenraeume plotten falls notwendig
    if norm(ERaum1) ~= 0
        plot(ERaum1(1,:),ERaum1(2,:), 'k','LineWidth',1.0)
    end
    if norm(ERaum2) ~= 0
        plot(ERaum2(1,:),ERaum2(2,:), 'k','LineWidth',1.0)
    end
    
    %Loesung plotten
    for j = 1:size(AWerte,2)
        [t,y] = ode45(fun4,[0 5*Farbe],[double(AWerte(1,j)) double(AWerte(2,j))]);
        plot(y(:,1),y(:,2), FC)
    end
    end
end


%%
%Funktionen

%!Doppelte reelle EW!
%Scherung, dh doppelter EW 0 (nicht diag)
function [A,B1,B2] = fkt_Scherung(VJ,xValue,yValue,k)
    
    EV1 = VJ(:,1)/norm(VJ(:,1));
    
    a1 = [xValue-k*EV1(1)+k*EV1(2); yValue-k*EV1(2)-EV1(1)];
    a2 = [xValue-2*k*EV1(1)+2*k*EV1(2); yValue-2*k*EV1(2)-2*k*EV1(1)];
    a3 = [xValue-3*k*EV1(1)+3*k*EV1(2); yValue-3*k*EV1(2)-3*k*EV1(1)];
    a4 = [xValue+k*EV1(1)-k*EV1(2); yValue+k*EV1(2)+EV1(1)];
    a5 = [xValue+2*k*EV1(1)-2*k*EV1(2); yValue+2*k*EV1(2)+2*k*EV1(1)];
    a6 = [xValue+3*k*EV1(1)-3*k*EV1(2); yValue+3*k*EV1(2)+3*k*EV1(1)];
    A = [a1 a2 a3 a4 a5 a6];

    B1 = [xValue-3*EV1(1) xValue+3*EV1(1); yValue-3*EV1(2) yValue+3*EV1(2)]; 
    B2 =  0;
    
end

%entartete Knoten, dh J nicht diag
function [A,B1,B2]=fkt_ent_Knoten(VJ,WJ,xValue,yValue,a,k,d)    
    EV1 = VJ(:,1)/norm(VJ(:,1));
    if WJ(1,1)>0
        a1 = [xValue+k*EV1(1)+d*EV1(2); yValue+k*EV1(2)-d*EV1(1)];
        a2 = [xValue+k*EV1(1)-d*EV1(2); yValue+k*EV1(2)+d*EV1(1)];
        a3 = [xValue-k*EV1(1)+d*EV1(2); yValue-k*EV1(2)-d*EV1(1)];
        a4 = [xValue-k*EV1(1)-d*EV1(2); yValue-k*EV1(2)+d*EV1(1)];
        a5 = [xValue+(a+k)*EV1(1)+d*EV1(2); yValue+(a+k)*EV1(2)-d*EV1(1)];
        a6 = [xValue+(a+k)*EV1(1)-d*EV1(2); yValue+(a+k)*EV1(2)+d*EV1(1)];
        a7 = [xValue-(a+k)*EV1(1)+d*EV1(2); yValue-(a+k)*EV1(2)-d*EV1(1)];
        a8 = [xValue-(a+k)*EV1(1)-d*EV1(2); yValue-(a+k)*EV1(2)+d*EV1(1)];
        A = [a1 a2 a3 a4 a5 a6 a7 a8];
    else
        c = 0.66*k;
        b = 0.66*a;
        a1 = [xValue+c*EV1(1)-d*EV1(2); yValue+c*EV1(2)+d*EV1(1)];
        a2 = [xValue+(b+c)*EV1(1)-d*EV1(2); yValue+(b+c)*EV1(2)+d*EV1(1)];
        a3 = [xValue+(2*b+c)*EV1(1)-d*EV1(2); yValue+(2*b+c)*EV1(2)+d*EV1(1)];
        a4 = [xValue-k*EV1(1)+d*EV1(2); yValue-c*EV1(2)-d*EV1(1)];
        a5 = [xValue-(b+c)*EV1(1)+d*EV1(2); yValue-(b+c)*EV1(2)-d*EV1(1)];
        a6 = [xValue-(2*b+c)*EV1(1)+d*EV1(2); yValue-(2*b+c)*EV1(2)-d*EV1(1)];
        A =  [a1 a2 a3 a4 a5 a6];
    end
    B1 = [xValue-3*EV1(1) xValue+3*EV1(1); yValue-3*EV1(2) yValue+3*EV1(2)]; 
    B2 =  0;
end

%Stern, dh J diag (zwei EV)
function [A, B1, B2] = fkt_Stern(VJ,xValue,yValue,d)
    EV1=(VJ(:,1)/norm(VJ(:,1)));
    EV2=(VJ(:,2)/norm(VJ(:,2)));
    
    a1 = [xValue+d*EV1(1)+d*EV2(1); yValue+d*EV1(2)+d*EV2(2)];
    a2 = [xValue+d*EV1(1)-d*EV2(1); yValue+d*EV1(2)-d*EV2(2)];
    a3 = [xValue-d*EV1(1)+d*EV2(1); yValue-d*EV1(2)+d*EV2(2)];
    a4 = [xValue-d*EV1(1)-d*EV2(1); yValue-d*EV1(2)-d*EV2(2)];

    A = [a1 a2 a3 a4];
    % Eigenräume
    B1 = [xValue-3*EV1(1) xValue+3*EV1(1); 
          yValue-3*EV1(2) yValue+3*EV1(2)]; 
    B2 = [xValue-3*EV2(1) xValue+3*EV2(1); 
          yValue-3*EV2(2) yValue+3*EV2(2)];
end

%!zwei verschiedene reele EW!
%Sattel
function [A,B1,B2] = fkt_Sattel(VJ,WJ,xValue,yValue,a,k,d)
    
    EV1=(VJ(:,1)/norm(VJ(:,1)));
    EV2=(VJ(:,2)/norm(VJ(:,2)));
    
    if WJ(1,1)>WJ(2,2)
        a1 = [xValue+k*EV1(1)+d*EV2(1); yValue+k*EV1(2)+d*EV2(2)];
        a2 = [xValue+k*EV1(1)-d*EV2(1); yValue+k*EV1(2)-d*EV2(2)];
        a3 = [xValue-k*EV1(1)+d*EV2(1); yValue-k*EV1(2)+d*EV2(2)];
        a4 = [xValue-k*EV1(1)-d*EV2(1); yValue-k*EV1(2)-d*EV2(2)];
        a5 = [xValue+(a+k)*EV1(1)+d*EV2(1); yValue+(a+k)*EV1(2)+d*EV2(2)];
        a6 = [xValue+(a+k)*EV1(1)-d*EV2(1); yValue+(a+k)*EV1(2)-d*EV2(2)];
        a7 = [xValue-(a+k)*EV1(1)+d*EV2(1); yValue-(a+k)*EV1(2)+d*EV2(2)];
        a8 = [xValue-(a+k)*EV1(1)-d*EV2(1); yValue-(a+k)*EV1(2)-d*EV2(2)];
    else
        a1 = [xValue+k*EV2(1)+d*EV1(1); yValue+k*EV2(2)+d*EV1(2)];
        a2 = [xValue+k*EV2(1)-d*EV1(1); yValue+k*EV2(2)-d*EV1(2)];
        a3 = [xValue-k*EV2(1)+d*EV1(1); yValue-k*EV2(2)+d*EV1(2)];
        a4 = [xValue-k*EV2(1)-d*EV1(1); yValue-k*EV2(2)-d*EV1(2)];
        a5 = [xValue+(a+k)*EV2(1)+d*EV1(1); yValue+(a+k)*EV2(2)+d*EV1(2)];
        a6 = [xValue+(a+k)*EV2(1)-d*EV1(1); yValue+(a+k)*EV2(2)-d*EV1(2)];
        a7 = [xValue-(a+k)*EV2(1)+d*EV1(1); yValue-(a+k)*EV2(2)+d*EV1(2)];
        a8 = [xValue-(a+k)*EV2(1)-d*EV1(1); yValue-(a+k)*EV2(2)-d*EV1(2)];
    end
    A = [a1 a2 a3 a4 a5 a6 a7 a8];
    B1 = [xValue-3*EV1(1) xValue+3*EV1(1); 
          yValue-3*EV1(2) yValue+3*EV1(2)]; 
    B2 = [xValue-3*EV2(1) xValue+3*EV2(1); 
          yValue-3*EV2(2) yValue+3*EV2(2)];
end

% Knoten
function [A,B1,B2]=fkt_Knoten(VJ,WJ,xValue,yValue,a,k,d)
    
    EV1=(VJ(:,1)/norm(VJ(:,1)));
    EV2=(VJ(:,2)/norm(VJ(:,2)));    

    if WJ(1,1)>WJ(2,2)
      a1 = [xValue+k*EV1(1)+d*EV2(1); yValue+k*EV1(2)+d*EV2(2)];
      a2 = [xValue+k*EV1(1)-d*EV2(1); yValue+k*EV1(2)-d*EV2(2)];
      a3 = [xValue-k*EV1(1)+d*EV2(1); yValue-k*EV1(2)+d*EV2(2)];
      a4 = [xValue-k*EV1(1)-d*EV2(1); yValue-k*EV1(2)-d*EV2(2)];
      a5 = [xValue+(a+k)*EV1(1)+d*EV2(1); yValue+(a+k)*EV1(2)+d*EV2(2)];
      a6 = [xValue+(a+k)*EV1(1)-d*EV2(1); yValue+(a+k)*EV1(2)-d*EV2(2)];
      a7 = [xValue-(a+k)*EV1(1)+d*EV2(1); yValue-(a+k)*EV1(2)+d*EV2(2)];
      a8 = [xValue-(a+k)*EV1(1)-d*EV2(1); yValue-(a+k)*EV1(2)-d*EV2(2)];
    else
      a1 = [xValue+k*EV2(1)+d*EV1(1); yValue+k*EV2(2)+d*EV1(2)];
      a2 = [xValue+k*EV2(1)-d*EV1(1); yValue+k*EV2(2)-d*EV1(2)];
      a3 = [xValue-k*EV2(1)+d*EV1(1); yValue-k*EV2(2)+d*EV1(2)];
      a4 = [xValue-k*EV2(1)-d*EV1(1); yValue-k*EV2(2)-d*EV1(2)];
      a5 = [xValue+(a+k)*EV2(1)+d*EV1(1); yValue+(a+k)*EV2(2)+d*EV1(2)];
      a6 = [xValue+(a+k)*EV2(1)-d*EV1(1); yValue+(a+k)*EV2(2)-d*EV1(2)];
      a7 = [xValue-(a+k)*EV2(1)+d*EV1(1); yValue-(a+k)*EV2(2)+d*EV1(2)];
      a8 = [xValue-(a+k)*EV2(1)-d*EV1(1); yValue-(a+k)*EV2(2)-d*EV1(2)];
    end
    A = [a1 a2 a3 a4 a5 a6 a7 a8];
    B1 = [xValue-3*EV1(1) xValue+3*EV1(1); 
          yValue-3*EV1(2) yValue+3*EV1(2)]; 
    B2 = [xValue-3*EV2(1) xValue+3*EV2(1); 
          yValue-3*EV2(2) yValue+3*EV2(2)];
end

%lineare Kontraktion/Expansion d.h. ein EW=0
function [A, B1, B2] = fkt_lin_Knoten(VJ1,VJ2, xWert, yWert, s, k)
    
    EV1 = VJ1/norm(VJ1);
    EV2 = VJ2/norm(VJ2);
    
    B1 = [xWert+3*EV1(1) xWert-3*EV1(1); yWert+3*EV1(2) yWert-3*EV1(2)];
    B2 = [xWert+3*EV2(1) xWert-3*EV2(1); yWert+3*EV2(2) yWert-3*EV2(2)];
    
        
    a1 = [xWert+s*EV1(1)+k*EV2(1); yWert+s*EV1(2)+k*EV2(2)];
    a2 = [xWert+s*EV1(1)-k*EV2(1); yWert+s*EV1(2)-k*EV2(2)];
    a3 = [xWert-s*EV1(1)+k*EV2(1); yWert-s*EV1(2)+k*EV2(2)];
    a4 = [xWert-s*EV1(1)-k*EV2(1); yWert-s*EV1(2)-k*EV2(2)];
    a5 = [xWert+2*s*EV1(1)+k*EV2(1); yWert+2*s*EV1(2)+k*EV2(2)];
    a6 = [xWert+2*s*EV1(1)-k*EV2(1); yWert+2*s*EV1(2)-k*EV2(2)];
    a7 = [xWert-2*s*EV1(1)+k*EV2(1); yWert-2*s*EV1(2)+k*EV2(2)];
    a8 = [xWert-2*s*EV1(1)-k*EV2(1); yWert-2*s*EV1(2)-k*EV2(2)];
    
    A = [a1 a2 a3 a4 a5 a6 a7 a8];
end

%!komplexe EW!
%Zentrum
function [A] = fkt_Zentrum(xValue,yValue,d)
    a1 = [xValue; yValue+d];
    a2 = [xValue; yValue+2*d];
    a3 = [xValue; yValue+3*d];
    a4 = [xValue; yValue+4*d];
    A = [a1 a2 a3 a4];
end

%Strudel
function [A] = fkt_Strudel(xValue,yValue,r)
        a1 = [xValue+r; yValue+r];
        a2 = [xValue-r; yValue-r];
        a3 = [xValue+r; yValue-r];
        a4 = [xValue-r; yValue+r];
        A = [a1 a2 a3 a4];
end


%%
%Klassifikation

function [A, B1, B2, R] = fkt_Klassifikation(J, xWert, yWert)
        
    [VJ,DJ] = eig(J);
    lambda = DJ(1,1);
    mu = DJ(2,2);
    
    %überprüfe, ob Jacobische ~=0
    if norm(J)==0
        A = 0;
        B1 = 0;
        B2 = 0;
        R = 0;    
    
    elseif mu == lambda
        if mu * lambda == 0
            A = fkt_Scherung(VJ, xWert, yWert, 0.66);
            B1 = 0;
            B2 = 0;
            R = 1;
            Sc=1
        else
            if lambda > 0
                if rank(VJ) == 1
                    %entartet, instab Knoten, nicht diagonalisierbar
                    [A, B1, B2] = fkt_ent_Knoten(VJ,DJ,xWert,yWert,1.5,1.5,2);
                    R = -1;
                    ieK
                else
                    %entartet, instab Knoten, diagonalisierbar 
                    [A, B1, B2] = fkt_Stern(VJ,xWert, yWert, 2);
                    R = -1;
                    iSt=1
                end
            else
                if rank(VJ) == 1
                    %entartet, stab Knoten, nicht diagonalisierbar
                    [A, B1, B2] = fkt_ent_Knoten(VJ,DJ,xWert,yWert,1.5,1.5,2);
                    R = 1;
                    seK=1
                else
                    %entartet, stab Knoten, diagonalisierbar
                    [A, B1, B2] = fkt_Stern(VJ,xWert, yWert, 2);
                    R = 1;
                    sSt=1
                end
            end
        end
     
    else %lambda ungleich mu
        if lambda * mu == 0
            if mu == 0
                if lambda > 0
                    %lineare Expansion mit mu = 0
                    [A, B1, B2] = fkt_lin_Knoten(VJ(:,2),VJ(:,1),xWert, yWert,1,2);
                    R = -1;
                else
                    %lineare Kontraktion
                    [A, B1, B2] = fkt_lin_Knoten(VJ(:,2),VJ(:,1),xWert, yWert,1,2);
                    R = 1;
                end
            else
                if mu > 0
                    %lineare Expansion mit lambda=0
                    [A, B1, B2] = fkt_lin_Konten(VJ(:,1),VJ(:,2),xWert, yWert,1,2);
                    R = -1;
                else 
                    %lineare Kontraktion
                    [A, B1, B2] = fkt_lin_Knoten(VJ(:,1),VJ(:,2),xWert, yWert,1,2);
                    R = 1;
                end
            end
            
        else
            if abs(imag(lambda)) < 0.001
                if lambda * mu < 0
                    [A, B1, B2] = fkt_Sattel(VJ,DJ,xWert,yWert,0.3,0.1,4);
                    R = 1.0001;
                    Sa=1
                else
                    if lambda > 0
                        %instabil
                        [A, B1, B2] = fkt_Knoten(VJ,DJ,xWert,yWert,1.5,1.5,2);
                        R = -1;
                        iK=1
                    else
                        %stabil
                        [A, B1, B2] = fkt_Knoten(VJ,DJ,xWert,yWert,1.5,1.5,2);
                        R = 1;
                        sK=1
                    end
                end
                
            else %lambda imaginaer
                if abs(real(lambda)) < 0.001
                    if lambda>0
                        A = fkt_Zentrum(xWert,yWert,2);
                        R = -1;
                        B1 = 0;
                        B2 = 0;
                    else
                        A = fkt_Zentrum(xWert,yWert,0.5);
                        R = 1;
                        B1 = 0;
                        B2 = 0;
                    end
                else
                    if real(lambda) > 0
                        %instabil
                        A = fkt_Strudel(xWert,yWert,2);
                        R = -1;
                        B1 = 0;
                        B2 = 0;
                        iS=1
                    else
                        %stabil
                        A = fkt_Strudel(xWert,yWert,2);
                        R = 1;
                        B1 = 0;
                        B2 = 0;
                        sS=1
                    end
                end
            end
        end
    end

end
%%
%Video
