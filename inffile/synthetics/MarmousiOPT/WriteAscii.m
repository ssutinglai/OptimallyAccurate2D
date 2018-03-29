clear all
nt=5500;
Size=nt*40; % due to the bug in OPT-FD
time=zeros(nt);
dt=10^(-8);

time=(1:nt)*10^(-8);


for iReceiver=1:39;
    iReceiver
    num5digit=sprintf('%05d',100+iReceiver*20);
    num2digit=sprintf('%02d',iReceiver);
    filename=strcat('00501.',num5digit,'.OPT_UX');
    filename2=strcat('Receiver',num2digit,'.txt');
    
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


for iReceiver=1:39;
    iReceiver
    num5digit=sprintf('%05d',100+iReceiver*20);
    num2digit=sprintf('%02d',iReceiver);
    filename=strcat('00501.',num5digit,'.OPT_UZ');
    filename2=strcat('Receiver',num2digit,'Z.txt');
    
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
