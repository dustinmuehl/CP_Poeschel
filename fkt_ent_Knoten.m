function [A,B1,B2]=fkt_ent_Knoten(VJ,WJ,xValue,yValue,k,r)    
    
    if WJ(1,1)<0
    c = 0.66*k;
    a1 = [xValue+c*VJ(1,1)+r yValue+c*VJ(2,1)+r];
    a2 = [xValue+2*c*VJ(1,1)+r yValue+2*c*VJ(2,1)+r];
    a3 = [xValue+3*c*VJ(1,1)+r yValue+3*c*VJ(2,1)+r];
    a4 = [xValue-c*VJ(1,1)-r yValue-c*VJ(2,1)-r];
    a5 = [xValue-2*c*VJ(1,1)-r yValue-2*c*VJ(2,1)-r];
    a6 = [xValue-3*c*VJ(1,1)-r yValue-3*c*VJ(2,1)-r];
    A = [a1; a2; a3; a4; a5; a6];
    else
    a1 = [xValue+k*VJ(1,1)+r yValue+k*VJ(2,1)+r];
    a2 = [xValue+k*VJ(1,1)-r yValue+k*VJ(2,1)-r];
    a3 = [xValue-k*VJ(1,1)+r yValue-k*VJ(2,1)+r];
    a4 = [xValue-k*VJ(1,1)-r yValue-k*VJ(2,1)-r];
    a5 = [xValue+2*k*VJ(1,1)+r yValue+2*k*VJ(2,1)+r];
    a6 = [xValue+2*k*VJ(1,1)-r yValue+2*k*VJ(2,1)-r];
    a7 = [xValue-2*k*VJ(1,1)+r yValue-2*k*VJ(2,1)+r];
    a8 = [xValue-2*k*VJ(1,1)-r yValue-2*k*VJ(2,1)-r];
    A = [a1; a2; a3; a4; a5; a6; a7; a8];
    end
    B1 = [xValue-3*VJ(1,1) xValue+3*VJ(1,1); 
          yValue-3*VJ(2,1) yValue+3*VJ(2,1)];
    B2 = 0;
end
