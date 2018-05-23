clear all;clc;close all


load ../../../../Katayama/iva1296/tek00281.csv
load ../../../../Katayama/iva1296/tek00282.csv

load ./Receiver01P.txt
load ./Receiver01S.txt
load ./Receiver01X.txt
load ./Receiver01Z.txt

load ../../waveletOptcopy

wt=waveletOptcopy(1:5000,1)*(2e-9)-(1.6e-6);
wu=waveletOptcopy(1:5000,2);

t277=tek00281(:,1);
u277=tek00281(:,2);
t278=tek00282(:,1);
u278=tek00282(:,2);

tp=Receiver01P(:,1);
up=Receiver01P(:,2);
ts=Receiver01S(:,1);
us=Receiver01S(:,2);

Xt=Receiver01X(:,1);
Xu=Receiver01X(:,2);
Zt=Receiver01Z(:,1);
Zu=Receiver01Z(:,2);

figure(1)
subplot(2,1,1)

plot(tp-1.6e-6,-up*0.08-5e-3,'b'); grid on; hold on
plot(t277,u277,'k')
axis([-inf 12*10^-6 -inf inf])
title('Strain synthetics for P-wave component')

subplot(2,1,2)
%plot(ts02,-us02-(-up02*0.944),'b'); grid on; hold on
plot(ts-1.6e-6,us,'b'); grid on; hold on
plot(t278,u278,'k')
axis([-inf 12*10^-6 -inf inf])
title('Strain synthetics for S-wave component')


% figure(2)
% subplot(2,1,1)
% 
% plot(Xt,Xu,'b'); grid on; hold on
% plot(t277,u277,'k')
% axis([2*10^-6 12*10^-6 -inf inf])
% title('X component')
% 
% subplot(2,1,2)
% plot(Zt,Zu,'b'); grid on; hold on
% plot(t278,u278,'k')
% axis([2*10^-6 12*10^-6 -inf inf])
% title('Z component')
