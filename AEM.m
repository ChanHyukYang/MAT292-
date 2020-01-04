function [x,y] = AEM(t0,tN,y0,h,deriv)
vals=t0:h:tN;
tval=zeros(1,length(vals));
yval=zeros(1,length(vals));
yval(1)=y0;
tol=1e-8;
tval(1)=t0;
counter=1;
while tval(counter)<tN
    m1=deriv((tval(counter)),yval(counter));
    Y=yval(counter)+(h)*(m1);
    Z1=yval(counter)+(h/2)*(m1);
    mZ=deriv(tval(counter)+(h/2),Z1);
    Z=Z1+(h/2)*(mZ);
    D=Z-Y;
    if (abs(D))<tol
        yval(counter+1)=Z+D;
        tval(counter+1)=h+tval(counter);
        counter=counter+1;
    else
        h = 0.9*h*min(max(tol/abs(D),0.3),2);
    end
end
y = yval;
x = tval(1:counter);