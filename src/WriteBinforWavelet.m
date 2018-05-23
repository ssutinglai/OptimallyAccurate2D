clear all;clc;close all
load Wavelet.txt
time=Wavelet(:,1);
disp=Wavelet(1:5:end,2);

[nt,nx] = size(disp);

fileID = fopen('wavelet.bin','w');
fwrite(fileID,disp,'double');
fclose(fileID);

fileID = fopen('wavelet.bin','r');
input1 = fread(fileID,[nt,nx],'double');
fclose(fileID)

figure
    plot(input1)