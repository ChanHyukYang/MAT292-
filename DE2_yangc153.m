function [t,y] = DE2_yangc153(t0,tN,y0,y1,h,p,q,g)    
    t = t0:h:tN;
    n = length(t);
    y = zeros(size(t));    
    y(1) = y0;
    y(2) = y(1) + h*y1;
    for i = 2:n-1
        y(i+1) = (g(t(i)) - q(t(i))*(y(i)) - p(t(i))*((y(i)) - y(i-1))/h) * h^2 + 2*y(i) - y(i-1);
    end
end