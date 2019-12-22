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