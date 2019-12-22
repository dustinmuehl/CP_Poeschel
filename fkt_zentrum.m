function [A] = fkt_zentrum(xValue,yValue,a)
    a1 = [xValue+a yValue];
    a2 = [xValue+2*a yValue];
    a3 = [xValue+3*a yValue];
    a4 = [xValue+4*a yValue];
    a5 = [xValue+5*a yValue];
    A = [a1; a2; a3; a4; a5];
end