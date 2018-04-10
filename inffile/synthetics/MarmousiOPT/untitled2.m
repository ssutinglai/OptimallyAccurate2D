clear all
load Receiver01P.txt
load Receiver01S.txt
load ../../waveletOpt
load waveletOpt
%load Receiver002P.txt
%load Receiver002S.txt
tp01=Receiver01P(:,1);
up01=Receiver01P(:,2);
ts01=Receiver01S(:,1);
us01=Receiver01S(:,2);
tw=waveletOpt(:,1)*4*10^(-9)-1.6*10^(-6);
uw=waveletOpt(:,2);
%tp02=Receiver002P(:,1);
%up02=Receiver002P(:,2);
%ts02=Receiver002S(:,1);
%us02=Receiver002S(:,2);

figure(1)
subplot(2,1,1)
plot(tp01,-up01)
axis([-inf inf -inf inf])
title('Strain synthetics for p-wave')

subplot(2,1,2)
% plot(tp02,-up02)
% axis([1*10^(-5) 14*10^(-6) -inf inf])

% figure(2)
%subplot(2,1,2)
plot(ts01,us01)
axis([-inf inf -inf inf])

% subplot(2,1,2)
% plot(ts02,us02)
% axis([1*10^(-5) 2*10^(-5) -inf inf])

figure(2)
plot(tw,uw)