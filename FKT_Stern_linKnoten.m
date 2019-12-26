%%
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

%%
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