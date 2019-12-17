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