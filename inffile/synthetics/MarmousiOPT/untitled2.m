clear all;clc;close all


load ../../../../Katayama/iva1296/tek00277.csv
load ../../../../Katayama/iva1296/tek00278.csv

load ../Con(1.3)/Receiver01P.txt
load ../Con(1.3)/Receiver01S.txt
load ../Con(1.3)/Receiver01X.txt
load ../Con(1.3)/Receiver01Z.txt

load ../../waveletOptcopy

wt=waveletOptcopy(1:5000,1)*(2e-9)-(1.6e-6);
wu=waveletOptcopy(1:5000,2);

tsk01=tek00277(:,1);
usk01=tek00277(:,2);
ts01=tek00278(:,1);
us01=tek00278(:,2);

tp02=Receiver01P(:,1);
up02=Receiver01P(:,2);
ts02=Receiver01S(:,1);
us02=Receiver01S(:,2);

Xt=Receiver01X(:,1);
Xu=Receiver01X(:,2);
Zt=Receiver01Z(:,1);
Zu=Receiver01Z(:,2);

figure(1)
subplot(2,1,1)

plot(tp02,-up02*6,'b'); grid on; hold on
plot(tsk01,usk01*100,'k')
axis([2*10^-6 12*10^-6 -inf inf])
title('Strain synthetics for P-wave component')

subplot(2,1,2)
%plot(ts02,-us02-(-up02*0.944),'b'); grid on; hold on
plot(ts02,us02*2.2,'b'); grid on; hold on
plot(ts01,us01,'k')
axis([2*10^-6 12*10^-6 -inf inf])
title('Strain synthetics for S-wave component')


figure(2)
subplot(2,1,1)

plot(Xt,Xu*100000000000,'b'); grid on; hold on
plot(tsk01,usk01,'k')
axis([2*10^-6 12*10^-6 -inf inf])
title('X component')

subplot(2,1,2)
%plot(ts02,-us02-(-up02*0.944),'b'); grid on; hold on
plot(Zt,Zu,'b'); grid on; hold on
plot(ts01,us01,'k')
axis([2*10^-6 12*10^-6 -inf inf])
title('Z component')
