clear all
fileID1=fopen('2d_circle.vp','r');
AA=fread(fileID1,'single'); %For reading kind(1.e0)
seismogram=AA(1,:);
fclose(fileID1);
