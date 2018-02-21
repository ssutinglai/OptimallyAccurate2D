load TBz.dat
load BBz.dat
load LBx.dat
load RBx.dat
%a=TBz(1,:);
%b=BBz(1,:);
a=LBx(1,:);
b=RBx(1,:);
x=1:900
plot(x,a,'*',x,b,'.')