function [t,x] =IEMVector(t0,tN,x0,h,deriv)
vals=t0:h:tN;
tval=zeros(1,length(vals));
xval=zeros(length(x0),length(vals));
xvalup=zeros(length(x0),length(vals));
for i=1:length(x0)
    xval(i,1)=x0(i);
    xvalup(i,1)=x0(i);
end
for i=1:length(vals)
    tval(i)=(i-1)*h+t0;
end
for i=2:length(vals)
    m10=deriv(1,[xval(1,i-1),xval(2,i-1)]);
    m20=deriv(2,[xval(1,i-1),xval(2,i-1)]);
    xvalup(1,i)=xval(1,i-1)+(h*m10(1));
    xvalup(2,i)=xval(2,i-1)+(h*m20(2));
    m11=deriv(1,[xvalup(1,i),xvalup(2,i)]);
    m21=deriv(2,[xvalup(1,i),xvalup(2,i)]);
    xval(1,i)=xval(1,i-1)+(h/2)*(m10(1)+m11(1));
    xval(2,i)=xval(2,i-1)+(h/2)*(m20(2)+m21(2));
end
x(1,:) = xval(1,:);
x(2,:) = xval(2,:);
t = tval;

