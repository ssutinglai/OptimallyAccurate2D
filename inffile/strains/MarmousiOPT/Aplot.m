clear all
nx=600;
nz=600;

%dt=2*10^(-9);
%t0=-1.6*(10^(-6));
%time=t0+(1:nt)*2*10^(-9);

    fileID1=fopen('strainD00340.00301.00110.OPT_dat','r');
    AA=fread(fileID1,'single'); 
    %seismogram=AA(1:nt);
    strain=reshape(AA,[nx,nz]);
   
    fclose(fileID1);
    figure(1);
    image(strain*1.e10);
    
    fileID1=fopen('strainS00340.00301.00110.OPT_dat','r');
    AA=fread(fileID1,'single'); 
    %seismogram=AA(1:nt);
    strain=reshape(AA,[nx,nz]);
   
    fclose(fileID1);
    figure(2);
    image(strain*1.e10);
    
    
