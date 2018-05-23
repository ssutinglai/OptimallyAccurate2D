clear all
nt=5000;
Size=nt*40; % due to the bug in OPT-FD
time=zeros(nt);
dt=4*10^(-9);
t0=-1.6*(10^(-6));
time=t0+(1:nt)*4*10^(-9);



for iReceiver=1:5001;
    iReceiver
    num4digit=sprintf('%04d',-1+iReceiver*1);
    num2digit=sprintf('%04d',iReceiver+1);
    filename=strcat('strainD0',num4digit,'.00301.00110.OPT_dat');
    filename2=strcat('Receiver',num2digit,'P.txt');
    
    fileID1=fopen(filename,'r');
    AA=fread(fileID1,'single'); %For reading kind(1.e0)
    seismogram=AA(1:nt);
    
    %AA1=(1:Size)*10^(-8)
    %AA1=AA1'
    fclose(fileID1);
    A=[time',seismogram];
    
    %
    fileID = fopen(filename2,'w');
    fprintf(fileID,'%e %e\n',A');
    fclose(fileID);
    iReceiver
end

for iReceiver=1:5001;
    iReceiver
    num4digit=sprintf('%04d',-1+iReceiver*1);
    num2digit=sprintf('%04d',iReceiver+1);
    filename=strcat('strainS0',num4digit,'.00301.00110.OPT_dat');
    filename2=strcat('Receiver',num2digit,'S.txt');
    
    fileID1=fopen(filename,'r');
    AA=fread(fileID1,'single'); %For reading kind(1.e0)
    seismogram=AA(1:nt);
    
    %AA1=(1:Size)*10^(-8)
    %AA1=AA1'
    fclose(fileID1);
    A=[time',seismogram];
    
    %
    fileID = fopen(filename2,'w');
    fprintf(fileID,'%e %e\n',A');
    fclose(fileID);
    iReceiver
end

