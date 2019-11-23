for xValue = -2 : 0.5 : 2
    for yValue = -2 : 0.5 : 2
        c = xValue;
        d = yValue;
        H = dsolve('Dy = y-x','y(c)=d', 'x');
        dH = diff(H);
        y = matlabFunction(H);
        dy = matlabFunction(dH);
        for i = -2:0.5:2
            m = dy(xValue,yValue,i);
            y0 = y(xValue,yValue,i);
            syms z
            t = matlabFunction((y0-m*i)+m*z);
            xt = i-0.2:0.2:i+0.2;
            hold on
            if m == 0
                plot(xt,t())
            else
                plot(xt,t(xt))
            end
            hold off
        end
    end
end