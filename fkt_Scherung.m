function A = fkt_Scherung(xValue,yValue,k)
    a1 = [xValue-k yValue+k];
    a2 = [xValue-2*k yValue+2*k];
    a3 = [xValue-3*k yValue+3*k];
    a4 = [xValue+k yValue-k];
    a5 = [xValue+2*k yValue-2*k];
    a6 = [xValue+3*k yValue-3*k];
    A = [a1; a2; a3; a4; a5; a6];
end