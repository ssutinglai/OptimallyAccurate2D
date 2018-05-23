clear all, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code is created to make sure the correction of the velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ../test_VP_3900_VS_2520/Receiver01P.txt
pt=Receiver01P(:,1)-1.6e-6;
pu=(-Receiver01P(:,2)-5e-1)*0.5;

load ../../../../Katayama/iva1296/tek00315.csv

rt=tek00315(:,1);
ru=tek00315(:,2)*150;

plot(pt,pu,'r'); hold on
plot(rt,ru,'k');
axis([2e-6 12e-6 -inf inf])