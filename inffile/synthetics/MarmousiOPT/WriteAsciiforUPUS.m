clear all
nt=7000;
Size=nt*40; % due to the bug in OPT-FD
time=zeros(nt);
dt=2*10^(-9);
%t0=-1.6*(10^(-6));
t0=0;
time=t0+(1:nt)*2*10^(-9)-2*10^(-9);



for iReceiver=1;
    iReceiver
    num5digit=sprintf('%05d',100+iReceiver*10);
    num2digit=sprintf('%02d',iReceiver);
    filename=strcat('00301.',num5digit,'.OPT_UP');
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


for iReceiver=1;
    iReceiver
    num5digit=sprintf('%05d',100+iReceiver*10);
    num2digit=sprintf('%02d',iReceiver);
    filename=strcat('00301.',num5digit,'.OPT_US');
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

for iReceiver=1;
    iReceiver
    num5digit=sprintf('%05d',100+iReceiver*400);
    num2digit=sprintf('%02d',iReceiver);
    filename=strcat('00301.',num5digit,'.OPT_UX');
    filename2=strcat('Receiver',num2digit,'X.txt');
    
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

for iReceiver=1;
    iReceiver
    num5digit=sprintf('%05d',100+iReceiver*400);
    num2digit=sprintf('%02d',iReceiver);
    filename=strcat('00301.',num5digit,'.OPT_UZ');
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
