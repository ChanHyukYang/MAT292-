function [x,y] =IEM(t0,tN,y0,h,deriv)
vals=t0:h:tN;
tval=zeros(1,length(vals));
yval=zeros(1,length(vals));
yval(1)=y0;
yvalup=zeros(1,length(vals));
yvalup(1)=y0;
for i=1:length(vals)
    tval(i)=(i-1)*h+t0;
end
for i=2:length(vals)
    m1=deriv(tval(i-1),yval(i-1));
    yvalup(i)=yval(i-1)+(h*m1);
    m2=deriv(tval(i),yvalup(i));
    yval(i)=yval(i-1)+(h/2)*(m2+m1);
end
y = yval;
x = tval;
    
    
    
    
    