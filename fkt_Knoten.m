function [A,B1,B2]=fkt_Knoten(VJ,WJ,xValue,yValue,k,r)
    if WJ(1,1)>WJ(2,2)
    a1 = [xValue+k*VJ(1,1)+r yValue+k*VJ(2,1)+r];
    a2 = [xValue+k*VJ(1,1)-r yValue+k*VJ(2,1)-r];
    a3 = [xValue-k*VJ(1,1)+r yValue-k*VJ(2,1)+r];
    a4 = [xValue-k*VJ(1,1)-r yValue-k*VJ(2,1)-r];
    a5 = [xValue+2*k*VJ(1,1)+r yValue+2*k*VJ(2,1)+r];
    a6 = [xValue+2*k*VJ(1,1)-r yValue+2*k*VJ(2,1)-r];
    a7 = [xValue-2*k*VJ(1,1)+r yValue-2*k*VJ(2,1)+r];
    a8 = [xValue-2*k*VJ(1,1)-r yValue-2*k*VJ(2,1)-r];
    else
    a1 = [xValue+k*VJ(1,2)+r yValue+k*VJ(2,2)+r];
    a2 = [xValue+k*VJ(1,2)-r yValue+k*VJ(2,2)-r];
    a3 = [xValue-k*VJ(1,2)+r yValue-k*VJ(2,2)+r];
    a4 = [xValue-k*VJ(1,2)-r yValue-k*VJ(2,2)-r];
    a5 = [xValue+2*k*VJ(1,2)+r yValue+2*k*VJ(2,2)+r];
    a6 = [xValue+2*k*VJ(1,2)-r yValue+2*k*VJ(2,2)-r];
    a7 = [xValue-2*k*VJ(1,2)+r yValue-2*k*VJ(2,2)+r];
    a8 = [xValue-2*k*VJ(1,2)-r yValue-2*k*VJ(2,2)-r];
    end
    A = [a1; a2; a3; a4; a5; a6; a7; a8];
    B1 = [xValue-3*VJ(1,1) xValue+3*VJ(1,1); 
          yValue-3*VJ(2,1) yValue+3*VJ(2,1)]; 
    B2 =  [xValue-3*VJ(1,2) xValue+3*VJ(1,2); 
          yValue-3*VJ(2,2) yValue+3*VJ(2,2)]; 
    
